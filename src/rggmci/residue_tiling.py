"""Residue-level tiling for RG-GMCI pairs, from DIAMOND alignments against MIBiG proteins (optional, report-only).

Why: ClusterBlast's hit table names which reference genes each region matches, but not which residues of those
genes. When an assembly break cuts one giant modular PKS or NRPS gene, both pieces match the same reference gene,
so the whole-gene tiling test reads them as two copies of one machinery (`OVERLAPPING_PARALOG`) and the pair is
demoted with `ST-PARALOG_no_complementarity_proof`. Aligning the proteins shows which stretch of the reference each
piece covers. Pieces of one broken gene cover different stretches; two copies of the same machinery cover the same
ones.

How, for one pair and one MIBiG cluster both regions already match in KnownClusterBlast:
1. Every protein of each region is aligned against the cluster's proteins (DIAMOND blastp, >= 30% identity,
   e <= 1e-10, the settings of the gap-rescue screen).
2. Each side is placed on the cluster on its own, best bitscore first. Every residue of a region protein and every
   residue of a cluster protein is used once, and an alignment counts only if it adds >= RESIDUE_MIN_NEW residues on
   both. So a PKS module is placed where it matches best, not against every similar module of the reference.
3. `residue_overlap` is the share of cluster residues both sides were placed on, over the smaller side's total.
   Per cluster: COMPLEMENTARY when each side has >= RESIDUE_MIN_SIDE residues placed and the overlap is at most
   RESIDUE_MAX_OVERLAP; OVERLAPPING when each side has that many and the overlap is larger; THIN otherwise.
4. Per pair, over all shared clusters: COMPLEMENTARY_RESIDUES (some cluster complementary, none overlapping),
   OVERLAPPING_RESIDUES (the reverse), MIXED_RESIDUES (both), THIN_RESIDUES (neither), NO_SHARED_MIBIG_REFERENCE (no
   shared cluster is in the DIAMOND database). The detail fields describe the top reference: the shared cluster
   both sides match best (the higher of the two sides' lower summed bitscore), whatever its class, so the headline is
   never picked for its answer. MIXED is common for modular PKS against class relatives at about 50% identity, where
   one module cannot be told from another; it is an honest "cannot tell", not a weak rescue.

The result is evidence fields on each tested pair. It never changes a score or a confidence. Homology is
similarity, not product identity, and no contigs are joined.

Needs the DIAMOND binary (`diamond`) and a DIAMOND database of MIBiG proteins whose sequence ids start with the
MIBiG accession (`BGC0000001|1`, `BGC0000001.5|…`). Without either, the pairs are left untested and the status
says why. Standard library only, so the standalone rggmci package ships this file unchanged.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path
from typing import Any, Iterable

from .ziputil import regular_file_names

RESIDUE_MIN_IDENTITY = 30      # DIAMOND --id
RESIDUE_MAX_EVALUE = "1e-10"   # DIAMOND --evalue
RESIDUE_MIN_NEW = 60           # an alignment counts only if it adds this many unused residues on both proteins
RESIDUE_MIN_SIDE = 300         # residues each side must have placed on a cluster before overlap is judged
RESIDUE_MAX_OVERLAP = 0.10     # shared placed residues / smaller side, at or below which the sides are complementary
RESIDUE_SCOPES = ("st_paralog", "all")
_ACC_RE = re.compile(r"^(BGC\d{7})")
_FEATURE_RE = re.compile(r"^ {5}(\S+)\s+")
_QUAL_RE = re.compile(r'^ {21}/(\w+)(?:=(.*))?$')


def find_diamond(explicit: str | os.PathLike | None = None) -> str | None:
    """The DIAMOND binary: an explicit path, then $RGGMCI_DIAMOND, then `diamond` on PATH."""
    for cand in (explicit, os.environ.get("RGGMCI_DIAMOND")):
        if cand and Path(cand).is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return shutil.which("diamond")


def mibig_accession(ref: str) -> str:
    """`BGC0000996.5` or `BGC0000996|3` -> `BGC0000996`; anything else -> ''."""
    m = _ACC_RE.match(ref or "")
    return m.group(1) if m else ""


def _cds_proteins(gbk_text: str) -> list[tuple[str, str]]:
    """(name, translation) for every CDS with a translation in one GenBank record text."""
    out: list[tuple[str, str]] = []
    in_features, cur, quals, key = False, None, {}, None

    def flush() -> None:
        if cur == "CDS" and quals.get("translation"):
            name = quals.get("locus_tag") or quals.get("protein_id") or quals.get("gene") or f"cds{len(out) + 1}"
            out.append((name, quals["translation"]))

    for line in gbk_text.splitlines():
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if not in_features:
            continue
        if line.startswith("ORIGIN") or line.startswith("//"):
            flush()
            in_features, cur, quals, key = False, None, {}, None
            continue
        fm = _FEATURE_RE.match(line)
        if fm and not line.startswith(" " * 21):
            flush()
            cur, quals, key = fm.group(1), {}, None
            continue
        qm = _QUAL_RE.match(line)
        if qm:
            key = qm.group(1)
            quals[key] = (qm.group(2) or "").strip('"')
        elif key and line.startswith(" " * 21):
            quals[key] += line.strip().strip('"') if key == "translation" else " " + line.strip().strip('"')
    flush()
    return [(n, re.sub(r"[^A-Za-z*]", "", s).rstrip("*")) for n, s in out]


def _region_gbk_name(names: set[str], bgc: Any) -> str:
    name = getattr(bgc, "source_gbk", "") or ""
    if name in names:
        return name
    tail = f".region{int(getattr(bgc, 'region_number', 0) or 0):03d}.gbk"
    contig = getattr(bgc, "contig", "")
    hits = sorted(n for n in names if n.endswith(tail) and contig and contig in Path(n).name)
    return hits[0] if len(hits) == 1 else ""


def region_proteins(zip_path: str | Path, bgcs: Iterable[Any], wanted: set[str]) -> dict[str, list[tuple[str, str]]]:
    """Proteins of the wanted regions, read from each region's own GenBank file inside the antiSMASH ZIP."""
    out: dict[str, list[tuple[str, str]]] = {}
    with zipfile.ZipFile(zip_path) as zf:
        names = set(regular_file_names(zf))
        for b in bgcs:
            if b.bgc_id not in wanted:
                continue
            name = _region_gbk_name(names, b)
            if name:
                out[b.bgc_id] = _cds_proteins(zf.read(name).decode("utf-8", "replace"))
    return out


def _run_diamond(diamond: str, db: str, query: Path, out: Path, threads: int) -> list[tuple]:
    subprocess.run([diamond, "blastp", "-q", str(query), "-d", str(db), "-o", str(out), "--threads", str(threads),
                    "--quiet", "--sensitive", "--id", str(RESIDUE_MIN_IDENTITY), "--evalue", RESIDUE_MAX_EVALUE,
                    "--max-target-seqs", "0", "--outfmt", "6", "qseqid", "sseqid", "pident", "qstart", "qend",
                    "sstart", "send", "bitscore"], check=True, capture_output=True)
    rows = []
    with open(out, newline="") as fh:
        for q, s, p, qs, qe, ss, se, bs in csv.reader(fh, delimiter="\t"):
            rows.append((q, s, float(p), int(qs), int(qe), int(ss), int(se), float(bs)))
    return rows


def place_side(hsps: Iterable[tuple]) -> list[tuple[str, int, int, float, float]]:
    """One side's placement on one cluster: (subject, from, to, identity, bitscore) per counted alignment.

    `hsps` are (query, subject, identity, qstart, qend, sstart, send, bitscore). Best bitscore first, every residue
    of a query protein and of a subject protein used once, and an alignment counts only if it adds
    RESIDUE_MIN_NEW unused residues on both.
    """
    q_used: dict[str, set[int]] = {}
    s_used: dict[str, set[int]] = {}
    placed = []
    for q, s, p, qs, qe, ss, se, bs in sorted(hsps, key=lambda h: (-h[7], h[0], h[1], h[3], h[5])):
        qa, qb = sorted((qs, qe))
        sa, sb = sorted((ss, se))
        qu, su = q_used.setdefault(q, set()), s_used.setdefault(s, set())
        q_new = [i for i in range(qa, qb + 1) if i not in qu]
        s_new = [i for i in range(sa, sb + 1) if i not in su]
        if len(q_new) < RESIDUE_MIN_NEW or len(s_new) < RESIDUE_MIN_NEW:
            continue
        qu.update(q_new)
        su.update(s_new)
        placed.append((s, sa, sb, p, bs))
    return placed


def _residues(placed: list[tuple]) -> set[tuple[str, int]]:
    return {(x[0], i) for x in placed for i in range(x[1], x[2] + 1)}


def compare_sides(placed_a: list, placed_b: list) -> dict[str, Any]:
    """Overlap of two sides' placements on one cluster, and the per-cluster class."""
    ra, rb = _residues(placed_a), _residues(placed_b)
    shared = len(ra & rb)
    small = min(len(ra), len(rb))
    overlap = round(shared / small, 3) if small else 0.0
    if len(ra) >= RESIDUE_MIN_SIDE and len(rb) >= RESIDUE_MIN_SIDE:
        cls = "COMPLEMENTARY" if overlap <= RESIDUE_MAX_OVERLAP else "OVERLAPPING"
    else:
        cls = "THIN"

    def ident(placed: list) -> float | str:
        w = sum(x[2] - x[1] + 1 for x in placed)
        return round(sum((x[2] - x[1] + 1) * x[3] for x in placed) / w, 1) if w else ""

    return {"residues_a": len(ra), "residues_b": len(rb), "residue_overlap": overlap, "class": cls,
            "joint_bitscore": round(min(sum(x[4] for x in placed_a), sum(x[4] for x in placed_b)), 1),
            "identity_a": ident(placed_a), "identity_b": ident(placed_b)}


def _stretches(placed: list, limit: int = 6) -> str:
    top = sorted(placed, key=lambda x: -(x[2] - x[1]))[:limit]
    return "; ".join(f"{x[0]}:{x[1]}-{x[2]}@{x[3]:.1f}" for x in sorted(top))


def _verdict(classes: list[str]) -> str:
    comp, over = "COMPLEMENTARY" in classes, "OVERLAPPING" in classes
    if comp and over:
        return "MIXED_RESIDUES"
    if comp:
        return "COMPLEMENTARY_RESIDUES"
    if over:
        return "OVERLAPPING_RESIDUES"
    return "THIN_RESIDUES"


def add_residue_tiling(result: dict[str, Any], zip_path: str | Path, bgcs: list[Any], reference_map: dict[str, Any],
                       diamond_db: str | os.PathLike, *, diamond: str | os.PathLike | None = None,
                       scope: str = "st_paralog", threads: int = 4, work_dir: str | os.PathLike | None = None) -> None:
    """Add residue-tiling fields to the pairs in `result["ranked_pairs"]` (in place). Never changes confidence.

    scope "st_paralog" tests the pairs the whole-gene paralog gate demoted; "all" tests every ranked pair.
    Untested pairs get `residue_tiling` = NOT_TESTED. `result["residue_tiling_status"]` says what ran.
    """
    if scope not in RESIDUE_SCOPES:
        raise ValueError(f"scope must be one of {RESIDUE_SCOPES}")
    pairs = result.get("ranked_pairs") or []
    todo = [p for p in pairs if scope == "all" or "ST-PARALOG" in str(p.get("acceptance_gate", ""))]
    for p in pairs:
        p.update({"residue_tiling": "NOT_TESTED", "residue_top_reference": "", "residue_top_class": "",
                  "residue_complementary_references": "", "residue_overlapping_references": "", "residue_overlap": "",
                  "residues_a": "", "residues_b": "", "residue_identity_a": "", "residue_identity_b": "",
                  "residue_references_tested": "", "residue_stretches_a": "", "residue_stretches_b": ""})
    status = {"scope": scope, "pairs_tested": 0, "diamond_db": str(diamond_db),
              "thresholds": {"identity": RESIDUE_MIN_IDENTITY, "evalue": RESIDUE_MAX_EVALUE,
                             "min_new": RESIDUE_MIN_NEW, "min_side": RESIDUE_MIN_SIDE,
                             "max_overlap": RESIDUE_MAX_OVERLAP},
              "interpretation_guard": "Residue tiling is homology evidence; it never changes confidence, joins "
                                      "contigs or names a product."}
    result["residue_tiling_status"] = status
    binary = find_diamond(diamond)
    if not binary:
        status["state"] = "DIAMOND_NOT_FOUND"
        return
    if not Path(diamond_db).is_file():
        status["state"] = "DIAMOND_DB_NOT_FOUND"
        return
    if not todo:
        status["state"] = "NO_PAIRS_IN_SCOPE"
        return

    kcb: dict[str, set[str]] = {}
    for r in reference_map.get("reference_records") or []:
        d = r if isinstance(r, dict) else vars(r)
        if d.get("db_kind") == "knownclusterblast" and mibig_accession(d.get("ref", "")):
            kcb.setdefault(d.get("bgc_id"), set()).add(mibig_accession(d["ref"]))
    wanted = {x for p in todo for x in (p["bgc_a"], p["bgc_b"])}
    prots = region_proteins(zip_path, bgcs, wanted)

    tmp = tempfile.TemporaryDirectory(prefix="rggmci_residue_") if work_dir is None else None
    wd = Path(tmp.name if tmp else work_dir)
    wd.mkdir(parents=True, exist_ok=True)
    try:
        faa = wd / "region_proteins.faa"
        with open(faa, "w") as fh:
            for bid in sorted(prots):
                for i, (name, seq) in enumerate(prots[bid], 1):
                    if seq:
                        fh.write(f">{bid}|{i}|{re.sub(r'[^A-Za-z0-9_.-]', '_', name)}\n{seq}\n")
        rows = _run_diamond(binary, str(diamond_db), faa, wd / "region_proteins_vs_db.tsv", threads)
    finally:
        if tmp:
            tmp.cleanup()
    by_region_ref: dict[tuple[str, str], list[tuple]] = {}
    in_db: set[str] = set()
    for h in rows:
        acc = mibig_accession(h[1])
        in_db.add(acc)
        by_region_ref.setdefault((h[0].split("|", 1)[0], acc), []).append(h)

    for p in todo:
        a, b = p["bgc_a"], p["bgc_b"]
        shared = sorted((kcb.get(a, set()) & kcb.get(b, set())) & in_db)
        p["residue_references_tested"] = len(shared)
        if not shared:
            p["residue_tiling"] = "NO_SHARED_MIBIG_REFERENCE"
            continue
        per_ref = []
        for acc in shared:
            pa, pb = place_side(by_region_ref.get((a, acc), [])), place_side(by_region_ref.get((b, acc), []))
            per_ref.append((acc, compare_sides(pa, pb), pa, pb))
        classes = [c["class"] for _acc, c, _pa, _pb in per_ref]
        p["residue_tiling"] = _verdict(classes)
        p["residue_complementary_references"] = classes.count("COMPLEMENTARY")
        p["residue_overlapping_references"] = classes.count("OVERLAPPING")
        acc, c, pa, pb = min(per_ref, key=lambda x: (-x[1]["joint_bitscore"], x[0]))
        p.update({"residue_top_reference": acc, "residue_top_class": c["class"], "residue_overlap": c["residue_overlap"],
                  "residues_a": c["residues_a"], "residues_b": c["residues_b"],
                  "residue_identity_a": c["identity_a"], "residue_identity_b": c["identity_b"],
                  "residue_stretches_a": _stretches(pa), "residue_stretches_b": _stretches(pb)})
        status["pairs_tested"] += 1
    status["state"] = "RAN"
