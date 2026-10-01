"""Reference-guided completion: each edge region's MIBiG reference, searched across the whole genome (report-only).

An assembly break can leave part of a pathway on a contig antiSMASH never flagged, or cut one gene in two. RG-GMCI
pairs edge regions that share ClusterBlast references; this module asks the reference itself where its genes are.

For every edge region (antiSMASH `contig_edge="True"`) whose best KnownClusterBlast hit is in the MIBiG protein
database, every protein of that reference is aligned against every protein of the genome, not only proteins inside
antiSMASH regions (one alignment call per genome covers all such references). Then:
- the gene table (`search`): per reference gene, PRESENT_IN_CORE, MISSING_FOUND_CLEAR, MISSING_FOUND_AMBIGUOUS or
  MISSING_NOT_FOUND, with the rules of tools/gap_directed_rescue.py, which imports this code;
- split genes (`split_genes`): one reference gene in two pieces at facing contig ends, called CLEAR, WEAK,
  RIVAL_STRONGER or MODULAR_UNRESOLVED; a piece pair that splits against references of unrelated compound families,
  at least one of them the reference of a core holding neither piece, is relabelled RECURRENT_COMMON_GENE (a broken
  common gene, such as a regulator, not pathway evidence). Two cores that each hold one piece and see the same split
  through their own references are the expected picture of a real split, not recurrence;
- partner tests, per contig other than the core that carries a clear find or a split piece:
  position (the find sits within PARTNER_END_BP of a contig end, or the contig is under PARTNER_SHORT_BP), read depth
  (SPAdes coverage ratio within DEPTH_MIN..1/DEPTH_MIN of the core's; blank when a name carries no coverage),
  reciprocal (the find's aligned stretch, searched against all MIBiG proteins, scores >= RECIPROCAL_MIN of its best
  score against the reference's compound family) and contig support (>= 2 finds that pass position and reciprocal,
  or one that is a piece of a CLEAR split gene). Mobile-element reference genes (transposase, integrase, recombinase,
  insertion element) never count toward support: they sit in many clusters and in no pathway. A contig is an
  accepted partner when support passes and depth does not fail.
- each RG-GMCI ranked pair gets `completion_tier`, `ref_completion_partner` (yes when one region's contig is an
  accepted partner for the other region's reference) and `split_gene_links` (the number of distinct genome piece pairs
  called CLEAR whose pieces sit on the two regions' contigs: one broken gene counts once, however many references see
  it).

Tiers, on every row and in the receipt: FULL, NO_ALIGNER, NO_MIBIG_PROTEINS, NO_WHOLE_GENOME_GENBANK, OFF. Anything
but FULL means the completion did not run; missing evidence is not absence of a partner. With no database the tier is
NO_MIBIG_PROTEINS whatever aligner is installed, so a run without a database reads the same on every machine. A ZIP
without a whole-genome GenBank file (region files only, as in a single-region download) cannot be searched outside its
regions, so its tier is NO_WHOLE_GENOME_GENBANK; it is decided from the input before any aligner is probed, and the
pairs are kept as RG-GMCI scored them.

Inputs are never searched for on disk: the MIBiG protein database comes from an explicit path, $RGGMCI_MIBIG_DB, or
(engine only) the workspace asset registry. The aligner is DIAMOND (explicit path, $RGGMCI_DIAMOND, PATH), else BLAST+
blastp and makeblastdb on PATH with the same thresholds; tables are comparable, not identical, and the aligner is
recorded per row.

Nothing here changes a score or a confidence. Homology is similarity, not product identity; no contigs are joined.
Standard library only, so the standalone rggmci package ships this file unchanged.
"""
from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Callable, Iterable

from ._gbk_shim import parse_genbank_text
from .csv_safety import SafeWriter
from .residue_tiling import find_diamond, mibig_accession
from .ziputil import regular_file_names

_log = logging.getLogger(__name__)

# ── Gene table and split-gene check (moved from tools/gap_directed_rescue.py, which imports them) ─────────────────
MIN_ID, MIN_COV, CLEAR_ID, CLEAR_RATIO = 30.0, 50.0, 35.0, 1.2
# inside the core region a weaker match still counts when the alignment is strong: the gene sits in the cluster, so its
# position supports the homology (a real nucleoside cluster had a phosphatase at 28-30% to its MIBiG homolog, e 1e-14)
MIN_ID_CORE, MAX_EVALUE_CORE = 25.0, 1e-10
PARTNER_MIN, CONCENTRATION, SHORT_CONTIG = 2, 0.6, 30000
BIOSYNTHETIC_KINDS = {"biosynthetic", "biosynthetic-additional"}
TABLE_COLS = ["reference_gene", "name", "reference_product", "reference_len_aa", "reference_gene_kind", "modular_pks", "status",
              "best_identity_pct", "best_coverage_pct", "reciprocal_best", "best_protein", "best_locus", "best_len_aa",
              "best_contig",
              "best_region_identity", "rivals", "bitscore_ratio_to_second"]
PARTNER_COLS = ["partner_contig", "partner_identity", "contig_length", "clear_finds", "biosynthetic_clear_finds",
                "genes", "concentrated", "split_plausible"]
SPLIT_MIN_AA, SPLIT_MAX_PIECE_COV, SPLIT_MAX_OVERLAP, SPLIT_MIN_UNION, SPLIT_END_BP = 40, 80.0, 20, 50.0, 300
SPLIT_RIVAL_MARGIN = 10.0  # identity points both pieces must beat a whole-gene match elsewhere by, for a CLEAR call
SPLIT_COLS = ["reference_gene", "name", "reference_len_aa", "status", "split_call", "modular_reference_gene",
              "modular_piece", "reference_union_pct", "overlap_aa",
              "piece1_locus", "piece1_region_identity", "piece1_reference_range", "piece1_identity_pct",
              "piece1_open_end_to_contig_end_bp",
              "piece2_locus", "piece2_region_identity", "piece2_reference_range", "piece2_identity_pct",
              "piece2_open_end_to_contig_end_bp",
              "whole_gene_rival_locus", "whole_gene_rival_identity_pct"]


def region_of(p: dict, regions: list) -> dict | None:
    return next((r for r in regions if r["contig"] == p["contig"] and r["start"] < p["end"] and r["end"] > p["start"]),
                None)


def search(ref, prots, regions, hits, core):
    """-> rows (one per reference gene) and partner rows."""
    in_core = lambda p: p["contig"] == core["contig"] and p["start"] < core["end"] and p["end"] > core["start"]

    def kept(h) -> bool:
        if h["sseqid"] not in prots or float(h["qcovhsp"]) < MIN_COV:
            return False
        if float(h["pident"]) >= MIN_ID:
            return True
        ev = h.get("evalue")  # precomputed tables may carry no e-value: then the 30% floor holds everywhere
        return (ev not in (None, "") and float(h["pident"]) >= MIN_ID_CORE and float(ev) <= MAX_EVALUE_CORE
                and in_core(prots[h["sseqid"]]))

    by_q = {}
    for h in hits:
        if kept(h):
            by_q.setdefault(h["qseqid"], []).append(h)
    for q, values in by_q.items():
        unique = {}
        for h in sorted(values, key=lambda h: (-float(h["bitscore"]), -float(h["pident"]),
                                              -float(h["qcovhsp"]), h["sseqid"])):
            unique.setdefault(h["sseqid"], h)
        by_q[q] = list(unique.values())
    rows = []
    for g in ref:
        hs = by_q.get(g["id"], [])
        row = {"reference_gene": g["i"], "name": g["name"], "reference_product": g["product"],
               "reference_len_aa": len(g["aa"]), "reference_gene_kind": g["kind"], "modular_pks": g["ks"] > 0,
               "status": "MISSING_NOT_FOUND"}
        inside = [h for h in hs if in_core(prots[h["sseqid"]])]
        outside = [h for h in hs if not in_core(prots[h["sseqid"]])]
        best = None
        if inside:
            best, row["status"] = inside[0], "PRESENT_IN_CORE"
        elif outside:
            best = outside[0]
            second = outside[1] if len(outside) > 1 else None
            clear = float(best["pident"]) >= CLEAR_ID and (second is None or
                                                          float(best["bitscore"]) >= CLEAR_RATIO * float(second["bitscore"]))
            row["status"] = "MISSING_FOUND_CLEAR" if clear else "MISSING_FOUND_AMBIGUOUS"
            row["rivals"] = len(outside) - 1
            row["bitscore_ratio_to_second"] = round(float(best["bitscore"]) / float(second["bitscore"]), 2) if second else ""
        if best:
            p = prots[best["sseqid"]]
            reg = region_of(p, regions)
            row.update(best_identity_pct=round(float(best["pident"]), 1), best_coverage_pct=round(float(best["qcovhsp"])),
                       best_protein=best["sseqid"], best_locus=p["tag"], best_len_aa=len(p["aa"]), best_contig=p["contig"],
                       best_region_identity=reg["identity"] if reg else f"{p['shown']} (no antiSMASH region)")
        rows.append(row)
    # modular PKS genes are left out: their best whole-gene match follows module paralogy (use KS placement instead)
    # reciprocal best: of the reference genes whose best match is one genome protein, only the highest-identity one is
    # that protein's own match. The others are cross-hits (typically paralogous PKS modules) and get no ribbon.
    top_for = {}
    for r in rows:
        if r.get("best_protein"):
            b = top_for.get(r["best_protein"])
            if b is None or r["best_identity_pct"] > b["best_identity_pct"]:
                top_for[r["best_protein"]] = r
    for r in rows:
        if r.get("best_protein"):
            r["reciprocal_best"] = top_for[r["best_protein"]] is r
    clear = [r for r in rows if r["status"] == "MISSING_FOUND_CLEAR" and r["best_contig"] != core["contig"]
             and not r["modular_pks"] and r.get("reciprocal_best")]
    total_clear = sum(1 for r in rows if r["status"] == "MISSING_FOUND_CLEAR"
                      and not r["modular_pks"] and r.get("reciprocal_best"))
    by_c = Counter(r["best_contig"] for r in clear)
    top = [c for c, _ in by_c.most_common(2)]
    concentrated = total_clear and sum(by_c[c] for c in top) / total_clear >= CONCENTRATION
    core_edge = core["edge"] in ("True", "true", "Edge")
    partners = []
    for c, n in by_c.most_common():
        if n < PARTNER_MIN:
            continue
        genes = [r for r in clear if r["best_contig"] == c]
        p0 = prots[genes[0]["best_protein"]]
        regs = [region_of(prots[r["best_protein"]], regions) for r in genes]
        partners.append({"partner_contig": c, "partner_identity": genes[0]["best_region_identity"],
                         "contig_length": p0["contig_len"], "clear_finds": n,
                         "biosynthetic_clear_finds": sum(1 for r in genes if r["reference_gene_kind"] in BIOSYNTHETIC_KINDS),
                         "genes": "; ".join(f"{r['reference_gene']} {r['name']} ({r['best_identity_pct']}%)" for r in genes),
                         "concentrated": bool(concentrated and c in top),
                         "split_plausible": bool(core_edge or p0["contig_len"] < SHORT_CONTIG
                                                 or any(r and r["edge"] in ("True", "true") for r in regs))})
    return rows, partners


def _open_end_distance(p: dict, terminus: str) -> int:
    """bp from a protein's N- or C-terminal end to the nearer end of its contig."""
    fwd = p["strand"] >= 0
    x = (p["start"] if fwd else p["end"]) if terminus == "N" else (p["end"] if fwd else p["start"])
    return min(x, p["contig_len"] - x)


def split_genes(ref, prots, regions, hits) -> tuple[list[dict], str]:
    """Reference genes found as two pieces at facing contig ends -> (rows, check status). Never changes the gene table."""
    if any("qstart" not in h or "qend" not in h for h in hits):
        return [], "not run: the hit table has no reference coordinates (qstart, qend)"
    best = {}
    for h in hits:
        k = (h["qseqid"], h["sseqid"])
        if h["sseqid"] in prots and float(h["pident"]) >= MIN_ID and (k not in best or
                                                                      float(h["bitscore"]) > float(best[k]["bitscore"])):
            best[k] = h
    by_q = {}
    for (q, _), h in sorted(best.items()):
        by_q.setdefault(q, []).append(h)
    span = lambda h: int(h["qend"]) - int(h["qstart"]) + 1
    where = lambda p: (region_of(p, regions) or {}).get("identity") or f"{p['shown']} (no antiSMASH region)"
    out = []
    for i, g in enumerate(ref):
        n_aa, hs = len(g["aa"]), by_q.get(g["id"], [])
        pieces = [h for h in hs if span(h) >= SPLIT_MIN_AA and float(h["qcovhsp"]) < SPLIT_MAX_PIECE_COV]
        found = None
        for x in range(len(pieces)):
            for y in range(x + 1, len(pieces)):
                a, b = sorted((pieces[x], pieces[y]), key=lambda h: (int(h["qstart"]), int(h["qend"])))
                pa, pb = prots[a["sseqid"]], prots[b["sseqid"]]
                if pa["contig"] == pb["contig"]:
                    continue
                overlap = max(0, min(int(a["qend"]), int(b["qend"])) - max(int(a["qstart"]), int(b["qstart"])) + 1)
                union = span(a) + span(b) - overlap
                # a covers the reference's start, so its C-terminal end is open; b covers the end, so its N-terminal end is
                da, db = _open_end_distance(pa, "C"), _open_end_distance(pb, "N")
                if (overlap > SPLIT_MAX_OVERLAP or 100 * union / n_aa < SPLIT_MIN_UNION
                        or da > SPLIT_END_BP or db > SPLIT_END_BP):
                    continue
                score = float(a["bitscore"]) + float(b["bitscore"])
                if found is None or score > found[0]:
                    found = (score, a, b, da, db, overlap, union)
        if not found:
            continue
        _, a, b, da, db, overlap, union = found
        pa, pb = prots[a["sseqid"]], prots[b["sseqid"]]
        rival = max((h for h in hs if float(h["qcovhsp"]) >= SPLIT_MAX_PIECE_COV
                     and h["sseqid"] not in (a["sseqid"], b["sseqid"])),
                    key=lambda h: float(h["bitscore"]), default=None)
        low = min(float(a["pident"]), float(b["pident"]))
        if g.get("modular") or pa.get("modular") or pb.get("modular"):
            call = "MODULAR_UNRESOLVED"   # module paralogy; order these with stretch assignment or KS placement
        elif rival is not None and float(rival["pident"]) >= low:
            call = "RIVAL_STRONGER"       # a whole gene elsewhere matches as well: likely paralog fragments
        elif rival is not None and low - float(rival["pident"]) < SPLIT_RIVAL_MARGIN:
            call = "WEAK"
        else:
            call = "CLEAR"
        out.append({"reference_gene": g["i"], "name": g["name"], "reference_len_aa": n_aa,
                    "status": "SPLIT_ACROSS_CONTIG_ENDS", "split_call": call,
                    "modular_reference_gene": bool(g.get("modular")),
                    "modular_piece": bool(pa.get("modular") or pb.get("modular")),
                    "reference_union_pct": round(100 * union / n_aa),
                    "overlap_aa": overlap,
                    "piece1_locus": pa["tag"], "piece1_region_identity": where(pa),
                    "piece1_reference_range": f"{a['qstart']}-{a['qend']}", "piece1_identity_pct": round(float(a["pident"]), 1),
                    "piece1_open_end_to_contig_end_bp": da,
                    "piece2_locus": pb["tag"], "piece2_region_identity": where(pb),
                    "piece2_reference_range": f"{b['qstart']}-{b['qend']}", "piece2_identity_pct": round(float(b["pident"]), 1),
                    "piece2_open_end_to_contig_end_bp": db,
                    "whole_gene_rival_locus": prots[rival["sseqid"]]["tag"] if rival else "",
                    "whole_gene_rival_identity_pct": round(float(rival["pident"]), 1) if rival else "",
                    "_i": i, "_pieces": (a["sseqid"], b["sseqid"])})
    return out, "run"


# ── Partner checks and reference discovery (moved from tools/gap_directed_rescue.py, which imports them) ──────────
# partner checks: is a find on another contig part of this pathway, or a housekeeping gene or a paralog-family member?
PARTNER_WINDOW_BP, PARALOG_LIMIT, RECIPROCAL_MARGIN, NEIGHBOUR_GENES, HOUSEKEEPING_MIN = 10000, 5, 5.0, 4, 2
# the reciprocal test: a find fails when its best score against the reference's compound family is below this share of
# its best MIBiG score (PKS and NRPS modules match many clusters almost equally, so a score ratio, not a rank)
RECIPROCAL_MIN = 0.9
# a short, named list of primary-metabolism operon markers (lysine, arginine, ribosome, tRNA charging, transcription,
# ATP synthase, gyrase); two of these beside a lone find mark a housekeeping block. Extend with evidence, not by guess.
HOUSEKEEPING_PFAM = re.compile(r"^(DHDPS|DapB_[NC]|ArgJ|Arg_repressor(_C)?|Arginosuc_synth|Arginosuc_syn_C|"
                               r"Semialdhyde_dh[C1_]*|Ribosomal_\w+|tRNA-synt_\w+|RNA_pol_\w+|ATP-synt\w*|"
                               r"DNA_gyraseB\w*|DNA_topoisoIV\w*)$")
PARTNER_CHECK_COLS = ["partner_verdict", "reciprocal_best_mibig", "reciprocal_best_identity", "reciprocal_reference_identity",
                      "paralogs_in_genome", "adjacent_finds", "housekeeping_neighbours", "neighbour_pfams",
                      "reciprocal_family_ratio"]
# reference discovery when none is given: KnownClusterBlast rank 1, else the MIBiG cluster with the most core-region
# proteins at >= 35% identity over >= 50% of the protein (at least 2, one of them biosynthetic), at most 250 kb
DISCOVERY_MIN_ID, DISCOVERY_MIN_COV, DISCOVERY_MIN_PROTEINS, DISCOVERY_MAX_KB = 35.0, 50.0, 2, 250.0
# when a database holds no cluster lengths, cluster size is judged by protein count instead (about 1 gene per kb)
DISCOVERY_MAX_PROTEINS = 250


def _pfam_names(seqs: dict[str, str], pfam_hmm: str, cpus: int = 4, only: re.Pattern | None = None
                ) -> dict[str, list[str]]:
    """Pfam domain names per protein (pyhmmer hmmscan against a pressed Pfam-A.hmm); {} when pyhmmer is absent.

    With `only`, just the profiles whose names match it are scanned (the housekeeping check needs only
    HOUSEKEEPING_PFAM), and E-values are still computed over every profile in the file (hmmscan Z), so the same domains
    pass the E <= 1e-5 cut as in a full scan."""
    try:
        import pyhmmer
        from pyhmmer.easel import Alphabet, TextSequence
        from pyhmmer.plan7 import HMMFile
    except ImportError:
        return {}
    aa = Alphabet.amino()
    digital = [TextSequence(name=k.encode(), sequence=v).digitize(aa) for k, v in seqs.items()]
    out: dict[str, list[str]] = {}
    with HMMFile(pfam_hmm) as hf:
        profiles = hf.optimized_profiles() if hf.is_pressed() else list(hf)
        extra = {}
        if only is not None:
            profiles = list(profiles)
            extra["Z"] = len(profiles)
            name = lambda p: p.name.decode() if isinstance(p.name, bytes) else p.name
            profiles = [p for p in profiles if only.match(name(p))]
        for top in pyhmmer.hmmer.hmmscan(digital, profiles, cpus=cpus, E=1e-5, **extra):
            q = top.query.name if hasattr(top, "query") else top.query_name
            q = q.decode() if isinstance(q, bytes) else q
            out[q] = [(h.name.decode() if isinstance(h.name, bytes) else h.name) for h in top if h.included]
    return out


def search_mibig(query_fasta: str, db: str | os.PathLike, threads: int = 1, sensitivity: str | None = None,
                 evalue: str = "1e-5", max_target_seqs: int = 25) -> dict[str, Any]:
    """Query proteins against a MIBiG protein database: a DIAMOND .dmnd, or a database folder (its .dmnd with DIAMOND,
    its FASTA with BLAST+). Same result shape as mamey.diamond_align.search_db: {ok, hits, reason}."""
    db = str(db)
    aligner, binary = find_aligner()
    folder = Path(db) if Path(db).is_dir() else None
    dmnd = db if db.endswith(".dmnd") else (str(folder / (DB_DMND + ".dmnd")) if folder else "")
    faa = str(folder / DB_FAA) if folder else ""
    if aligner == "blastp" and not faa:
        return {"ok": False, "hits": [], "reason": "BLAST+ needs a database folder with its FASTA"}
    if not aligner:
        return {"ok": False, "hits": [], "reason": "no aligner (DIAMOND, or BLAST+ blastp + makeblastdb)"}
    if aligner == "diamond" and not Path(dmnd).is_file():
        return {"ok": False, "hits": [], "reason": f"no DIAMOND database: {dmnd}"}
    with tempfile.TemporaryDirectory(prefix="rggmci_search_") as td:
        try:
            hits = align(aligner, binary, Path(query_fasta), Path(faa or dmnd), Path(td), threads=threads,
                         sensitivity=sensitivity or "default", subject_db=dmnd if aligner == "diamond" else "",
                         evalue=evalue, max_target_seqs=max_target_seqs)
        except (subprocess.CalledProcessError, OSError) as exc:
            return {"ok": False, "hits": [], "reason": f"{aligner} search failed: {type(exc).__name__}"}
    return {"ok": True, "hits": hits, "reason": f"ok ({aligner})"}


def _partner_finds(rows, prots, hits, core, ref=None):
    """The finds partner_checks judges (outside the core), the sequence each is searched with (the stretch the reference
    aligned to when the hit has genome coordinates, else the whole protein), and the reference-gene -> query id map."""
    gid = {g["i"]: g["id"] for g in ref} if ref else {}
    qid = lambda i: gid.get(int(i), f"g{int(i):03d}")
    in_core = lambda p: p["contig"] == core["contig"] and p["start"] < core["end"] and p["end"] > core["start"]
    finds = [r for r in rows if r.get("status") in ("MISSING_FOUND_CLEAR", "MISSING_FOUND_AMBIGUOUS")
             and r.get("best_protein") in prots and not in_core(prots[r["best_protein"]])]
    span = {}
    for h in hits:
        if "sstart" in h and "send" in h:
            k = (h["qseqid"], h["sseqid"])
            if k not in span or float(h["bitscore"]) > span[k][0]:
                span[k] = (float(h["bitscore"]), *sorted((int(h["sstart"]), int(h["send"]))))
    seqs = {}
    for r in finds:
        s = span.get((qid(r["reference_gene"]), r["best_protein"]))
        aa = prots[r["best_protein"]]["aa"]
        seqs[r["best_protein"]] = aa[s[1] - 1:s[2]] if s else aa
    return finds, seqs, qid


def partner_checks(rows, prots, hits, core, ref_acc, mibig_db=None, pfam_hmm=None, threads=4, sensitivity=None, *,
                   search: Callable | None = None, pfam_names: Callable | None = None,
                   compounds: dict[str, list[str]] | None = None, ref: list[dict] | None = None) -> dict:
    """Evidence for every find outside the core: is it part of this pathway, a housekeeping gene, or a paralog?

    - reciprocal: the find's MIBiG matches, every target (needs mibig_db). It fails when its best score against the
      reference's compound family (same accession, or a compound sharing a name stem, from `compounds`) is under RECIPROCAL_MIN of
      its best MIBiG score. The aligned stretch is searched when the hit carries genome coordinates (sstart, send), else
      the whole protein. The best cluster's identity and the reference's own identity are reported beside the ratio.
    - paralogs: other genome proteins matching the same reference gene (>= 30% over >= 50%).
    - adjacent finds: other finds of this reference on the same contig within 10 kb.
    - housekeeping neighbours: primary-metabolism Pfam domains among the 4 genes either side (needs pfam_hmm + pyhmmer).
    Verdicts: PARALOG_FAMILY (reciprocal fails, or >= 5 paralogs), SUPPORTED (an adjacent find), HOUSEKEEPING_CONTEXT
    (a lone find with >= 2 housekeeping neighbours), SINGLE_GENE (a lone find otherwise). Writes the fields into rows and
    returns a summary. Nothing here changes a status; the figure draws SUPPORTED and SINGLE_GENE contigs only.
    `search` and `pfam_names` default to search_mibig and _pfam_names; `ref` (reference gene dicts) maps reference gene
    numbers to query ids, which otherwise follow the tool's g001 form."""
    search = search or search_mibig
    pfam_names = pfam_names or _pfam_names
    finds, seqs, qid = _partner_finds(rows, prots, hits, core, ref)
    per_q = {}
    for h in hits:
        if h["sseqid"] in prots and float(h["pident"]) >= MIN_ID and float(h["qcovhsp"]) >= MIN_COV:
            per_q.setdefault(h["qseqid"], set()).add(h["sseqid"])
    recip, recip_status = {}, "not run: no MIBiG DIAMOND database (--mibig-db)"
    if mibig_db and finds:
        with tempfile.TemporaryDirectory() as td:
            qf = Path(td, "finds.faa")
            qf.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
            # every MIBiG target, not a top-N: a conserved family (SARP regulators, ABC transporters) matches many
            # clusters almost equally, and a cut-off can hide the reference family's own hit and fail a real gene
            res = search(str(qf), str(mibig_db), threads=threads, sensitivity=sensitivity, max_target_seqs=0)
        if res.get("ok"):
            recip_status = "run"
            for h in res["hits"]:
                acc = re.match(r"(BGC\d{7})", h["sseqid"])
                if not acc:
                    continue
                d = recip.setdefault(h["qseqid"], {"best": None, "ref": None, "fam": 0.0})
                cand = (float(h["bitscore"]), acc.group(1), float(h["pident"]))
                if d["best"] is None or cand > d["best"]:
                    d["best"] = cand
                if acc.group(1) == ref_acc and (d["ref"] is None or cand > d["ref"]):
                    d["ref"] = cand
                if same_family(acc.group(1), ref_acc, compounds or {}):
                    d["fam"] = max(d["fam"], cand[0])
        else:
            recip_status = f"not run: {res.get('reason')}"
    ref_in_db = bool(re.match(r"BGC\d{7}$", ref_acc or ""))
    if recip_status == "run" and not ref_in_db:
        recip_status = "run; not judged: the reference is not a MIBiG entry"
    lone = []
    for r in finds:
        p = prots[r["best_protein"]]
        r["paralogs_in_genome"] = max(0, len(per_q.get(qid(r["reference_gene"]), set())) - 1)
        r["adjacent_finds"] = sum(1 for o in finds if o is not r and prots[o["best_protein"]]["contig"] == p["contig"]
                                  and abs(prots[o["best_protein"]]["start"] - p["start"]) <= PARTNER_WINDOW_BP)
        d = recip.get(r["best_protein"])
        if d and d["best"]:
            r["reciprocal_best_mibig"], r["reciprocal_best_identity"] = d["best"][1], round(d["best"][2], 1)
            r["reciprocal_reference_identity"] = round(d["ref"][2], 1) if d["ref"] else ""
            r["reciprocal_family_ratio"] = round(d["fam"] / d["best"][0], 3) if d["best"][0] else ""
        # the reciprocal test can only judge when the reference itself is a MIBiG entry in the database
        fails = bool(ref_in_db and d and d["best"] and d["fam"] < RECIPROCAL_MIN * d["best"][0])
        if fails or r["paralogs_in_genome"] >= PARALOG_LIMIT:
            r["partner_verdict"] = "PARALOG_FAMILY"
        elif r["adjacent_finds"]:
            r["partner_verdict"] = "SUPPORTED"
        else:
            r["partner_verdict"] = "SINGLE_GENE"
            lone.append(r)
    pfam_status = "not run: no Pfam-A.hmm (--pfam)" if lone else "not needed"
    if pfam_hmm and lone:
        neigh = neighbour_ids(prots, [r["best_protein"] for r in lone])
        wanted = {k: prots[k]["aa"] for ks in neigh.values() for k in ks}
        names = pfam_names(wanted, str(pfam_hmm), threads) if wanted else {}
        pfam_status = "run" if names or not wanted else "not run: pyhmmer unavailable"
        for r in lone:
            nb = neigh[r["best_protein"]]
            r["neighbour_pfams"] = "; ".join(f"{prots[k]['tag']}:{','.join(names.get(k, [])[:2]) or '-'}" for k in nb)
            r["housekeeping_neighbours"] = sum(1 for k in nb if any(HOUSEKEEPING_PFAM.match(n) for n in names.get(k, [])))
            if r["housekeeping_neighbours"] >= HOUSEKEEPING_MIN:
                r["partner_verdict"] = "HOUSEKEEPING_CONTEXT"
    return {"reciprocal_mibig": recip_status, "neighbour_pfam": pfam_status,
            "verdicts": dict(Counter(r["partner_verdict"] for r in finds))}


def neighbour_ids(prots: dict, pids: Iterable[str]) -> dict[str, list[str]]:
    """The NEIGHBOUR_GENES genome proteins either side of each protein on its own contig (the housekeeping check's
    window)."""
    order = sorted(prots, key=lambda k: (prots[k]["contig"], prots[k]["start"]))
    pos = {k: i for i, k in enumerate(order)}
    out = {}
    for pid in pids:
        i, c = pos[pid], prots[pid]["contig"]
        out[pid] = [order[j] for j in range(max(0, i - NEIGHBOUR_GENES), i + NEIGHBOUR_GENES + 1)
                    if j < len(order) and j != i and prots[order[j]]["contig"] == c]
    return out


def _locus_kb(gbk: Path) -> float:
    first = open(gbk, encoding="utf-8", errors="replace").readline()
    m = re.search(r"(\d+)\s+bp", first)
    return int(m.group(1)) / 1000 if m else 0.0


def discover_references(queries: dict, prots: dict, mibig_db, mibig_dir=None, threads=4, sensitivity=None, hits=None,
                        *, search: Callable | None = None, lengths_kb: dict[str, float] | None = None,
                        protein_counts: dict[str, int] | None = None) -> dict:
    """{key: [protein ids]} -> {key: {reference, accession, proteins, bitscore, runner_up, tied, size_check}} by one
    search of every query protein against the MIBiG database (or precomputed `hits` with qseqid "key|protein id").

    Cluster size comes from `lengths_kb`, else the GenBank LOCUS line in `mibig_dir`, else (size_check
    "protein_count") at most DISCOVERY_MAX_PROTEINS proteins from `protein_counts`. `reference` is the GenBank path
    when `mibig_dir` is given, else the accession."""
    alias = {k: f"k{i}" for i, k in enumerate(queries)}   # DIAMOND keeps only the first word of a query name
    back = {v: k for k, v in alias.items()}
    if hits is None:
        search = search or search_mibig
        with tempfile.TemporaryDirectory() as td:
            qf = Path(td, "regions.faa")
            qf.write_text("".join(f">{alias[k]}|{pid}\n{prots[pid]['aa']}\n" for k, pids in queries.items() for pid in pids))
            res = search(str(qf), str(mibig_db), threads=threads, sensitivity=sensitivity, max_target_seqs=50)
        if not res.get("ok"):
            raise RuntimeError(f"reference discovery failed: {res.get('reason')}")
        hits = res["hits"]
    found = {}
    for h in hits:
        if float(h["pident"]) < DISCOVERY_MIN_ID or float(h["qcovhsp"]) < DISCOVERY_MIN_COV:
            continue
        key, pid = h["qseqid"].split("|", 1)
        key = back.get(key, key)
        acc = re.match(r"(BGC\d{7})", h["sseqid"])
        if acc:
            d = found.setdefault(key, {}).setdefault(acc.group(1), {})
            d[pid] = max(d.get(pid, 0.0), float(h["bitscore"]))

    def size_ok(acc: str) -> tuple[bool, str] | None:
        if lengths_kb and acc in lengths_kb:
            return lengths_kb[acc] <= DISCOVERY_MAX_KB, "length"
        if mibig_dir:
            gbk = Path(mibig_dir) / f"{acc}.gbk"
            return (_locus_kb(gbk) <= DISCOVERY_MAX_KB, "length") if gbk.exists() else None
        if protein_counts and acc in protein_counts:
            return protein_counts[acc] <= DISCOVERY_MAX_PROTEINS, "protein_count"
        return None

    out = {}
    for key, accs in found.items():
        ok, how = {}, {}
        for acc, got in accs.items():
            s = size_ok(acc)
            if s is None:
                continue
            if (len(got) >= DISCOVERY_MIN_PROTEINS and s[0]
                    and any(prots[p].get("kind") in BIOSYNTHETIC_KINDS for p in got)):
                ok[acc], how[acc] = got, s[1]
        ranked = sorted(ok.items(), key=lambda kv: (-len(kv[1]), -sum(kv[1].values()), kv[0]))
        if ranked:
            acc, got = ranked[0]
            runner = ranked[1] if len(ranked) > 1 else None
            out[key] = {"reference": Path(mibig_dir) / f"{acc}.gbk" if mibig_dir else acc, "accession": acc,
                        "proteins": len(got), "bitscore": round(sum(got.values())),
                        "runner_up": f"{runner[0]} ({len(runner[1])})" if runner else "",
                        "tied": bool(runner and len(runner[1]) == len(got)), "size_check": how[acc]}
    return out


# ── Completion settings ─────────────────────────────────────────────────────────────────────────────────────────────
TIERS = ("FULL", "NO_ALIGNER", "NO_MIBIG_PROTEINS", "NO_WHOLE_GENOME_GENBANK", "OFF")
SENSITIVITIES = ("sensitive", "more-sensitive", "ultra-sensitive", "default")
DEFAULT_SENSITIVITY = "sensitive"   # measured: DIAMOND's fast default misses 25-40% homologues; --sensitive finds them
ALIGN_MAX_EVALUE = "1e-5"           # no identity floor at alignment time; the table rules above decide what counts
PARTNER_END_BP, PARTNER_SHORT_BP = 20000, 40000
DEPTH_MIN = 0.67
SUPPORT_MIN_FINDS = 2
_MOBILE_RE = re.compile(r"transposase|integrase|recombinase|resolvase|insertion (?:sequence|element)|\bIS\d|phage",
                        re.I)


def is_mobile_gene(g: dict) -> bool:
    """A mobile-element reference gene, from its product or name text."""
    return bool(_MOBILE_RE.search(f"{g.get('product', '')} {g.get('name', '')}"))
MIBIG_DB_ENV = "RGGMCI_MIBIG_DB"
MIBIG_DB_ASSET = "mibig_local2088_proteins"
PFAM_ENV = "RGGMCI_PFAM_HMM"
PFAM_ASSET = "pfam_a_hmm"
DB_FAA, DB_TSV, DB_CLUSTERS, DB_MANIFEST, DB_DMND = ("mibig_proteins.faa", "mibig_proteins.tsv", "mibig_clusters.tsv",
                                                     "MANIFEST.json", "mibig_proteins")
GENE_TSV_COLS = ["id", "accession", "index", "name", "product", "kind", "ks", "modular", "length_aa"]
CLUSTER_TSV_COLS = ["accession", "compounds", "definition", "proteins", "length_bp"]
_COMMON_COLS = ["completion_tier", "aligner", "core_identity", "reference", "reference_source", "reference_compounds"]
COMPLETION_GENE_COLS = _COMMON_COLS + TABLE_COLS + PARTNER_CHECK_COLS
COMPLETION_SPLIT_COLS = _COMMON_COLS + SPLIT_COLS
COMPLETION_PARTNER_COLS = _COMMON_COLS + [
                           "partner_contig", "partner_region_identity", "contig_length", "finds", "passing_finds",
                           "genes", "position_ok", "depth_ratio", "depth_ok", "partner_verdicts",
                           "support_ok", "accepted", "reason"]
_NODE_RE = re.compile(r"^(NODE_\d+)_")


def _node_key(name: str) -> str:
    m = _NODE_RE.match(name or "")
    return m.group(1) if m else (name or "")


def family_stem(name: str) -> str:
    """'abyssomicin c' -> 'abyssomicin', so neoabyssomicin and atrop-abyssomicin count as one family."""
    words = re.findall(r"[a-z]{5,}", (name or "").lower())
    return max(words, key=len) if words else ""


def same_family(acc_a: str, acc_b: str, compounds: dict[str, list[str]]) -> bool:
    if acc_a == acc_b:
        return True
    stems_a = {family_stem(c) for c in compounds.get(acc_a, [])} - {""}
    stems_b = {family_stem(c) for c in compounds.get(acc_b, [])} - {""}
    return any(x in y or y in x for x in stems_a for y in stems_b)


# ── The MIBiG protein database ─────────────────────────────────────────────────────────────────────────────────────
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _gene_kind(q: dict) -> str:
    if q.get("gene_kind"):
        return q["gene_kind"][0]
    gf = " ".join(q.get("gene_functions", []))
    return re.split(r"[ (]", gf)[0] if gf else "other"


def reference_genes_from_gbk_text(text: str, accession: str) -> tuple[list[dict], str]:
    """Reference proteins of one MIBiG GenBank record, in coordinate order, with their assembly-line flag."""
    return _reference_record(text, accession)[:2]


def _reference_record(text: str, accession: str) -> tuple[list[dict], str, int]:
    """(genes, DEFINITION, record length in bp) of one MIBiG GenBank record."""
    recs = parse_genbank_text(text)
    if not recs:
        return [], "", 0
    rec = recs[0]
    cds = sorted([f for f in rec.features if f.type == "CDS" and f.qualifiers.get("translation")],
                 key=lambda f: f.location.start)
    doms = [(f.location.start, f.qualifiers.get("aSDomain", [""])[0]) for f in rec.features if f.type == "aSDomain"]
    genes = []
    for i, c in enumerate(cds, 1):
        q = c.qualifiers
        mine = [d for s, d in doms if c.location.start <= s < c.location.end]
        n_ks = mine.count("PKS_KS")
        modular = n_ks > 0 or any(d.startswith("Condensation") for d in mine) or mine.count("AMP-binding") >= 2
        genes.append({"id": f"{accession}|{i}", "accession": accession, "i": i,
                      "name": (q.get("gene") or q.get("locus_tag") or q.get("protein_id") or [f"g{i}"])[0],
                      "product": (q.get("product") or [""])[0], "kind": _gene_kind(q), "ks": n_ks, "modular": modular,
                      "aa": re.sub(r"[^A-Za-z]", "", q["translation"][0])})
    return genes, rec.description, len(rec.seq)


def build_mibig_db(gbk_dir: str | Path, out_dir: str | Path, *, compounds_json: str | Path | None = None,
                   diamond: str | os.PathLike | None = None) -> dict[str, Any]:
    """Write a MIBiG protein database from a folder of MIBiG GenBank files (one cluster per file).

    Files: mibig_proteins.faa (ids `<accession>|<n>`), mibig_proteins.tsv (per protein), mibig_clusters.tsv (compound
    names: from `compounds_json`, a list of {accession, compounds} under "entries", else none), MANIFEST.json, and
    mibig_proteins.dmnd when DIAMOND is available. Refuses a non-empty output folder. Never downloads anything.
    """
    src, out = Path(gbk_dir), Path(out_dir)
    files = sorted(p for p in src.glob("*.gbk") if not p.name.startswith("._"))
    if not files:
        raise FileNotFoundError(f"no .gbk files in {src}")
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"refusing to write into a non-empty folder: {out}")
    out.mkdir(parents=True, exist_ok=True)
    compounds: dict[str, list[str]] = {}
    if compounds_json:
        data = json.loads(Path(compounds_json).read_text(encoding="utf-8"))
        for e in data.get("entries", data if isinstance(data, list) else []):
            if e.get("accession"):
                compounds[mibig_accession(e["accession"]) or e["accession"]] = list(e.get("compounds") or [])
    n_prot = 0
    with open(out / DB_FAA, "w") as faa, open(out / DB_TSV, "w", newline="") as gt, \
            open(out / DB_CLUSTERS, "w", newline="") as ct:
        gw, cw = SafeWriter(gt, delimiter="\t"), SafeWriter(ct, delimiter="\t")
        gw.writerow(GENE_TSV_COLS)
        cw.writerow(CLUSTER_TSV_COLS)
        for f in files:
            acc = mibig_accession(f.stem) or f.stem
            genes, definition, length_bp = _reference_record(f.read_text(encoding="utf-8", errors="replace"), acc)
            for g in genes:
                faa.write(f">{g['id']}\n{g['aa']}\n")
                gw.writerow([g["id"], acc, g["i"], g["name"], g["product"], g["kind"], g["ks"], int(g["modular"]),
                             len(g["aa"])])
            n_prot += len(genes)
            cw.writerow([acc, "; ".join(compounds.get(acc, [])), definition, len(genes), length_bp])
    binary = find_diamond(diamond)
    dmnd = ""
    if binary:
        subprocess.run([binary, "makedb", "--in", str(out / DB_FAA), "-d", str(out / DB_DMND), "--quiet"],
                       check=True, capture_output=True)
        dmnd = DB_DMND + ".dmnd"
    manifest = {"tool": "rggmci build-mibig-db", "source_folder": src.name, "clusters": len(files), "proteins": n_prot,
                "compounds_source": Path(compounds_json).name if compounds_json else "",
                "files": {n: _sha256(out / n) for n in (DB_FAA, DB_TSV, DB_CLUSTERS)},
                "diamond_db": dmnd,
                "claim_safety": "MIBiG similarity is homology evidence, not product identity."}
    (out / DB_MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest


def load_mibig_db(path: str | Path) -> dict[str, Any]:
    """A MIBiG protein database folder -> genes by accession, compound names, and the file paths."""
    root = Path(path)
    if not (root / DB_FAA).is_file() or not (root / DB_TSV).is_file():
        raise FileNotFoundError(f"not a MIBiG protein database (needs {DB_FAA} and {DB_TSV}): {root}")
    seqs: dict[str, str] = {}
    name = None
    for line in open(root / DB_FAA, encoding="utf-8"):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    genes: dict[str, list[dict]] = {}
    with open(root / DB_TSV, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.setdefault(r["accession"], []).append(
                {"id": r["id"], "accession": r["accession"], "i": int(r["index"]), "name": r["name"],
                 "product": r["product"], "kind": r["kind"], "ks": int(r["ks"]), "modular": r["modular"] == "1",
                 "aa": seqs.get(r["id"], "")})
    compounds: dict[str, list[str]] = {}
    lengths_kb: dict[str, float] = {}
    if (root / DB_CLUSTERS).is_file():
        with open(root / DB_CLUSTERS, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                compounds[r["accession"]] = [c.strip() for c in r["compounds"].split(";") if c.strip()]
                if (r.get("length_bp") or "").isdigit():      # databases built before v9.7.445's length column: none
                    lengths_kb[r["accession"]] = int(r["length_bp"]) / 1000
    dmnd = root / (DB_DMND + ".dmnd")
    manifest = json.loads((root / DB_MANIFEST).read_text()) if (root / DB_MANIFEST).is_file() else {}
    return {"root": str(root), "faa": str(root / DB_FAA), "dmnd": str(dmnd) if dmnd.is_file() else "",
            "genes": genes, "compounds": compounds, "lengths_kb": lengths_kb,
            "protein_counts": {a: len(g) for a, g in genes.items()}, "manifest": manifest}


def _registry_path(asset: str, registry_from: str | os.PathLike) -> tuple[Path | None, str]:
    """The path of `asset` in the first OFFICIAL_DATA/ASSET_REGISTRY.tsv above `registry_from`."""
    for d in [Path(registry_from).resolve(), *Path(registry_from).resolve().parents]:
        reg = d / "OFFICIAL_DATA" / "ASSET_REGISTRY.tsv"
        if reg.is_file():
            for line in reg.read_text(encoding="utf-8", errors="replace").splitlines():
                p = line.split("\t")
                if len(p) >= 3 and p[0] == asset:
                    return d / p[2], "registry"
            return None, "not_registered"
    return None, "no_registry"


def resolve_pfam_hmm(explicit: str | os.PathLike | None = None, *, registry_from: str | os.PathLike | None = None
                     ) -> tuple[str | None, str]:
    """(Pfam-A.hmm path or None, how it was found) for the housekeeping-neighbour check. Order: explicit path,
    $RGGMCI_PFAM_HMM, then the asset registry row PFAM_ASSET above `registry_from` (engine only). Never searches the
    disk. The check also needs pyhmmer; without it the per-core status says so."""
    for how, cand in (("explicit", explicit), ("env", os.environ.get(PFAM_ENV))):
        if cand:
            return (str(cand), how) if Path(cand).is_file() else (None, f"{how}_path_missing")
    if registry_from:
        path, how = _registry_path(PFAM_ASSET, registry_from)
        if path is not None:
            return (str(path), "registry") if path.is_file() else (None, "registry_path_missing")
        return None, how
    return None, "not_given"


def resolve_mibig_db(explicit: str | os.PathLike | None = None, *, registry_from: str | os.PathLike | None = None
                     ) -> tuple[str | None, str]:
    """(database folder or None, how it was found). Order: explicit path, $RGGMCI_MIBIG_DB, then the asset registry
    (`OFFICIAL_DATA/ASSET_REGISTRY.tsv`, row MIBIG_DB_ASSET) of the first workspace above `registry_from`, when given.
    Never searches the disk for a database."""
    for how, cand in (("explicit", explicit), ("env", os.environ.get(MIBIG_DB_ENV))):
        if cand:
            return (str(cand), how) if Path(cand, DB_FAA).is_file() else (None, f"{how}_path_is_not_a_database")
    if registry_from:
        for d in [Path(registry_from).resolve(), *Path(registry_from).resolve().parents]:
            reg = d / "OFFICIAL_DATA" / "ASSET_REGISTRY.tsv"
            if reg.is_file():
                for line in reg.read_text(encoding="utf-8", errors="replace").splitlines():
                    p = line.split("\t")
                    if len(p) >= 3 and p[0] == MIBIG_DB_ASSET:
                        cand = d / p[2]
                        return (str(cand), "registry") if (cand / DB_FAA).is_file() else (None, "registry_path_missing")
                return None, "not_registered"
    return None, "not_given"


# ── Aligners ────────────────────────────────────────────────────────────────────────────────────────────────────────
_OUTFMT = ["qseqid", "sseqid", "pident", "qcovhsp", "bitscore", "qstart", "qend", "evalue", "sstart", "send"]


def find_aligner(diamond: str | os.PathLike | None = None) -> tuple[str, str]:
    """('diamond', path), ('blastp', path) when blastp and makeblastdb are both on PATH, or ('', '')."""
    d = find_diamond(diamond)
    if d:
        return "diamond", d
    b, m = shutil.which("blastp"), shutil.which("makeblastdb")
    return ("blastp", b) if b and m else ("", "")


def _write_faa(path: Path, items: Iterable[tuple[str, str]]) -> int:
    n = 0
    with open(path, "w") as fh:
        for k, s in items:
            if s:
                fh.write(f">{k}\n{s}\n")
                n += 1
    return n


def align(aligner: str, binary: str, query_faa: Path, subject_faa: Path, work: Path, *, threads: int = 4,
          sensitivity: str = DEFAULT_SENSITIVITY, subject_db: str = "", evalue: str = ALIGN_MAX_EVALUE,
          max_target_seqs: int = 0) -> list[dict]:
    """Every HSP at e <= `evalue` (all targets unless `max_target_seqs`), no identity floor, as dicts with query and
    subject coordinates."""
    out = work / f"hits_{aligner}_{query_faa.stem}.tsv"
    if aligner == "diamond":
        db = subject_db[:-5] if subject_db.endswith(".dmnd") else ""
        if not db:
            db = str(work / f"db_{subject_faa.stem}")
            subprocess.run([binary, "makedb", "--in", str(subject_faa), "-d", db, "--quiet"], check=True,
                           capture_output=True)
        mode = [] if sensitivity == "default" else [f"--{sensitivity}"]
        subprocess.run([binary, "blastp", "-q", str(query_faa), "-d", db, "-o", str(out), "--threads", str(threads),
                        "--quiet", *mode, "--evalue", str(evalue), "--max-target-seqs", str(max_target_seqs),
                        "--outfmt", "6", *_OUTFMT], check=True, capture_output=True)
    elif aligner == "blastp":
        mk = shutil.which("makeblastdb") or "makeblastdb"
        db = str(work / f"bdb_{subject_faa.stem}")
        subprocess.run([mk, "-in", str(subject_faa), "-dbtype", "prot", "-out", db], check=True, capture_output=True)
        subprocess.run([binary, "-query", str(query_faa), "-db", db, "-out", str(out), "-evalue", str(evalue),
                        "-max_target_seqs", str(max_target_seqs or 100000), "-num_threads", str(threads),
                        "-outfmt", "6 " + " ".join(_OUTFMT)], check=True, capture_output=True)
    else:
        raise ValueError(f"unknown aligner {aligner!r}")
    rows = []
    with open(out, newline="") as fh:
        for p in csv.reader(fh, delimiter="\t"):
            rows.append({"qseqid": p[0], "sseqid": p[1], "pident": float(p[2]), "qcovhsp": float(p[3]),
                         "bitscore": float(p[4]), "qstart": int(p[5]), "qend": int(p[6]), "evalue": float(p[7]),
                         "sstart": int(p[8]), "send": int(p[9])})
    return rows


# ── The genome ──────────────────────────────────────────────────────────────────────────────────────────────────────
def whole_genome_gbk_name(zf: zipfile.ZipFile) -> str:
    """The whole-genome GenBank file in an antiSMASH ZIP (the largest .gbk that is not a region file), or ''."""
    gbks = [n for n in regular_file_names(zf) if n.endswith(".gbk") and ".region" not in Path(n).name]
    return max(gbks, key=lambda n: zf.getinfo(n).file_size) if gbks else ""


def has_whole_genome_gbk(zip_path: str | Path) -> bool:
    try:
        with zipfile.ZipFile(zip_path) as zf:
            return bool(whole_genome_gbk_name(zf))
    except (OSError, zipfile.BadZipFile) as exc:   # an unreadable input is reported, not searched
        _log.debug("reference completion: cannot read %s as a ZIP: %s", zip_path, exc)
        return False


def genome_from_zip(zip_path: str | Path, label: str, bgcs: list[Any], edge_ids: set[str]) -> tuple[dict, list]:
    """Every CDS of the whole-genome GenBank in the antiSMASH ZIP (proteins), and the antiSMASH regions.

    Contigs are named as the region records name them (the bound source); a contig without a region takes its name
    from the ZIP's input FASTA, else from the GenBank record. SPAdes names are matched on their NODE number, because a
    record name can be shortened; other names must match exactly.
    """
    with zipfile.ZipFile(zip_path) as zf:
        names = regular_file_names(zf)
        genome = whole_genome_gbk_name(zf)
        if not genome:   # complete() records NO_WHOLE_GENOME_GENBANK before calling this; a direct caller still learns why
            raise FileNotFoundError("no whole-genome GenBank file in the ZIP")
        full: dict[str, str] = {}
        for fa in [n for n in names if n.endswith((".fasta", ".fa", ".fna")) and "input" in n]:
            for line in zf.read(fa).decode("utf-8", "replace").splitlines():
                if line.startswith(">"):
                    nm = line[1:].split()[0]
                    full.setdefault(_node_key(nm), nm)
        text = zf.read(genome).decode("utf-8", "replace")
    for b in bgcs:
        full[_node_key(b.contig)] = b.contig
    regions = [{"contig": b.contig, "start": int(b.start or 0), "end": int(b.end or 0),
                "n": int(b.region_number or 0), "edge": "True" if b.bgc_id in edge_ids else "False",
                "bgc_id": b.bgc_id,
                "identity": f"{label} / {b.contig} / region{int(b.region_number or 0):03d} / {b.bgc_id}"} for b in bgcs]
    prots: dict[str, dict] = {}
    for rec in parse_genbank_text(text):
        contig = full.get(_node_key(rec.name or rec.id), full.get(_node_key(rec.id), rec.name or rec.id))
        for f in rec.features:
            if f.type != "CDS" or not f.qualifiers.get("translation"):
                continue
            pid = f"q{len(prots) + 1:06d}"
            doms = [d.split(" (")[0] for d in f.qualifiers.get("sec_met_domain", [])]
            prots[pid] = {"contig": contig, "shown": f"{label} / {contig}", "start": f.location.start,
                          "end": f.location.end, "strand": f.location.strand or 1,
                          "tag": (f.qualifiers.get("locus_tag") or [pid])[0], "kind": _gene_kind(f.qualifiers),
                          "aa": re.sub(r"[^A-Za-z]", "", f.qualifiers["translation"][0]), "contig_len": len(rec.seq),
                          "modular": any(d == "PKS_KS" or d.startswith("Condensation") for d in doms)
                          or doms.count("AMP-binding") >= 2}
    return prots, regions


# ── Partner tests, recurrent split genes, and the join onto RG-GMCI pairs ──────────────────────────────────────────
def partner_tests(core: dict, ref: list[dict], rows: list[dict], splits: list[dict], prots: dict, regions: list,
                  hits: list[dict], depth_of: Callable[[str, str], float | None]) -> list[dict]:
    """One row per contig other than the core's carrying a find (MISSING_FOUND_CLEAR or _AMBIGUOUS, with its
    partner_checks verdict already on the row) or a split piece.

    A find passes when its reference gene is not a mobile element, it sits within PARTNER_END_BP of a contig end (or
    the contig is under PARTNER_SHORT_BP), and its verdict is SUPPORTED or SINGLE_GENE (PARALOG_FAMILY and
    HOUSEKEEPING_CONTEXT never pass). A CLEAR split piece passes on position alone. Support: >= SUPPORT_MIN_FINDS
    passing finds, one of them SUPPORTED, or a passing CLEAR split piece. Accepted: support, and depth does not fail."""
    by_gene = {g["i"]: g for g in ref}
    best_hit = {}
    for h in hits:
        k = (h["qseqid"], h["sseqid"])
        if k not in best_hit or h["bitscore"] > best_hit[k]["bitscore"]:
            best_hit[k] = h
    finds: dict[str, list[dict]] = {}
    for r in rows:
        if (r["status"] in ("MISSING_FOUND_CLEAR", "MISSING_FOUND_AMBIGUOUS") and r.get("best_contig")
                and r["best_contig"] != core["contig"]):
            g = by_gene[r["reference_gene"]]
            finds.setdefault(r["best_contig"], []).append({"gene": g, "protein": r["best_protein"],
                                                           "hit": best_hit.get((g["id"], r["best_protein"])),
                                                           "verdict": r.get("partner_verdict", "")})
    for s in splits:
        g = by_gene[s["reference_gene"]]
        for pid in s["_pieces"]:
            if prots[pid]["contig"] != core["contig"]:
                finds.setdefault(prots[pid]["contig"], []).append(
                    {"gene": g, "protein": pid, "hit": best_hit.get((g["id"], pid)), "split": s["split_call"]})
    out = []
    for contig, fs in finds.items():
        seen, uniq = set(), []
        for f in fs:
            k = (f["gene"]["id"], f["protein"], f.get("split", ""))
            if k not in seen:
                seen.add(k)
                uniq.append(f)
        p0 = prots[uniq[0]["protein"]]
        clen = p0["contig_len"]
        near_end = lambda p: min(p["start"], clen - p["end"]) <= PARTNER_END_BP or clen < PARTNER_SHORT_BP
        passing = [f for f in uniq if not is_mobile_gene(f["gene"]) and near_end(prots[f["protein"]])
                   and (f.get("split") == "CLEAR" or f.get("verdict") in ("SUPPORTED", "SINGLE_GENE"))]
        found = [f for f in passing if not f.get("split")]
        split_piece = any(f.get("split") == "CLEAR" for f in passing)
        support = (len(found) >= SUPPORT_MIN_FINDS and any(f["verdict"] == "SUPPORTED" for f in found)) or split_piece
        dr = depth_of(core["contig"], contig)
        dok = None if dr is None else dr >= DEPTH_MIN
        accepted = support and dok is not False
        reason = "accepted" if accepted else ("depth" if support else "no passing finds" if not passing else "support")
        reg = region_of(p0, regions)
        verdicts = Counter(f.get("verdict") or f"split {f['split']}" for f in uniq)
        out.append({"partner_contig": contig,
                    "partner_region_identity": reg["identity"] if reg else f"{p0['shown']} (no antiSMASH region)",
                    "contig_length": clen, "finds": len(uniq), "passing_finds": len(passing),
                    "genes": "; ".join(f"{f['gene']['i']} {f['gene']['name']}"
                                       + (f" ({f['hit']['pident']:.1f}%)" if f["hit"] else "")
                                       + (f" [{f['verdict']}]" if f.get("verdict") else "")
                                       + (f" [split {f['split']}]" if f.get("split") else "")
                                       + (" [mobile element]" if is_mobile_gene(f["gene"]) else "") for f in uniq),
                    "position_ok": any(near_end(prots[f["protein"]]) for f in uniq),
                    "depth_ratio": "" if dr is None else round(dr, 3), "depth_ok": "" if dok is None else dok,
                    "partner_verdicts": "; ".join(f"{k}:{v}" for k, v in sorted(verdicts.items())),
                    "support_ok": support, "accepted": accepted, "reason": reason})
    return sorted(out, key=lambda r: (not r["accepted"], -r["passing_finds"], r["partner_contig"]))


def mark_recurrent_splits(splits: list[dict], compounds: dict[str, list[str]]) -> int:
    """Relabel RECURRENT_COMMON_GENE a piece pair that splits against references of unrelated compound families when at
    least one of those references belongs to a core holding neither piece. `splits` rows carry `reference`,
    `_core_bgc`, `_pieces` and `_piece_bgcs` (the regions holding the pieces). Returns how many rows were relabelled."""
    seen: dict[frozenset, set[tuple[str, str]]] = {}
    for s in splits:
        seen.setdefault(frozenset(s["_pieces"]), set()).add((s["reference"], s["_core_bgc"]))
    n = 0
    for s in splits:
        both = sorted(seen[frozenset(s["_pieces"])])
        refs = sorted({r for r, _ in both})
        unrelated = any(not same_family(x, y, compounds) for i, x in enumerate(refs) for y in refs[i + 1:])
        foreign = any(core not in s["_piece_bgcs"] for _, core in both)
        if unrelated and foreign and s["split_call"] != "RECURRENT_COMMON_GENE":
            s["split_call_before_recurrence"] = s["split_call"]
            s["split_call"] = "RECURRENT_COMMON_GENE"
            n += 1
    return n


def join_pairs(pairs: list[dict], tier: str, cores: dict[str, dict], splits: list[dict], partners: list[dict],
               prots: dict) -> None:
    """Add completion_tier, ref_completion_partner and split_gene_links to every ranked pair (in place)."""
    accepted = {(p["_core_bgc"], p["partner_contig"]) for p in partners if p["accepted"]}
    for pr in pairs:
        pr["completion_tier"] = tier
        if tier != "FULL":
            pr["ref_completion_partner"] = ""
            pr["split_gene_links"] = ""
            continue
        a, b = pr.get("bgc_a"), pr.get("bgc_b")
        ca, cb = pr.get("contig_a"), pr.get("contig_b")
        pr["ref_completion_partner"] = "yes" if ((a, cb) in accepted or (b, ca) in accepted) else (
            "no" if (a in cores or b in cores) else "no_reference")
        links = set()
        for s in splits:
            if s["split_call"] != "CLEAR" or s["_core_bgc"] not in (a, b):
                continue
            if {prots[x]["contig"] for x in s["_pieces"]} == {ca, cb}:
                links.add(frozenset(s["_pieces"]))   # one broken gene, however many references see it
        pr["split_gene_links"] = len(links)


def _core_reference(reference_map: dict, edge_ids: set[str], db_genes: dict) -> dict[str, str]:
    """Each edge region's best-ranked KnownClusterBlast MIBiG accession that the database holds."""
    best: dict[str, tuple[int, str]] = {}
    for r in reference_map.get("reference_records") or []:
        d = r if isinstance(r, dict) else vars(r)
        acc = mibig_accession(d.get("ref", ""))
        if d.get("db_kind") != "knownclusterblast" or d.get("bgc_id") not in edge_ids or acc not in db_genes:
            continue
        rank = int(d.get("rank") or 10 ** 6)
        if d["bgc_id"] not in best or (rank, acc) < best[d["bgc_id"]]:
            best[d["bgc_id"]] = (rank, acc)
    return {k: v[1] for k, v in sorted(best.items())}


def complete(result: dict[str, Any], zip_path: str | Path, bgcs: list[Any], reference_map: dict[str, Any],
             edge_ids: set[str], *, label: str | None = None, mibig_db: str | os.PathLike | None = None,
             mibig_db_how: str = "", mode: str = "auto", diamond: str | os.PathLike | None = None,
             sensitivity: str = DEFAULT_SENSITIVITY, threads: int = 4,
             depth_of: Callable[[str, str], float | None] | None = None,
             work_dir: str | os.PathLike | None = None, pfam_hmm: str | os.PathLike | None = None,
             pfam_found_by: str = "") -> dict[str, Any]:
    """Run reference-guided completion and attach it to `result` (in place): result["reference_completion"] and three
    fields on every ranked pair. mode "off" records tier OFF and runs nothing. Returns the completion block."""
    if sensitivity not in SENSITIVITIES:
        raise ValueError(f"sensitivity must be one of {SENSITIVITIES}")
    label = label or Path(zip_path).stem
    if mode not in ("auto", "off"):
        raise ValueError("mode must be 'auto' or 'off'")
    have_db = mode != "off" and bool(mibig_db) and Path(mibig_db, DB_FAA).is_file()
    # the input is judged before the machine: a ZIP of region files only reads the same everywhere (see Tiers)
    has_genome = have_db and has_whole_genome_gbk(zip_path)
    aligner, binary = find_aligner(diamond) if has_genome else ("", "")   # probed only when it could be used
    if mode == "off":
        tier = "OFF"
    elif not have_db:
        tier = "NO_MIBIG_PROTEINS"
    elif not has_genome:
        tier = "NO_WHOLE_GENOME_GENBANK"
    else:
        tier = "FULL" if aligner else "NO_ALIGNER"
    block: dict[str, Any] = {"completion_tier": tier, "aligner": aligner, "sensitivity": sensitivity if aligner else "",
                             "mibig_db_found_by": mibig_db_how, "pfam_found_by": pfam_found_by,
                             "neighbour_pfam_scope": "housekeeping-marker profiles only (HOUSEKEEPING_PFAM); E-values over "
                                                     "every profile in the file", "cores": [], "genes": [], "splits": [], "partners": [],
                             "thresholds": {"align_max_evalue": ALIGN_MAX_EVALUE, "min_id": MIN_ID, "min_cov": MIN_COV,
                                            "clear_id": CLEAR_ID, "core_min_id": MIN_ID_CORE,
                                            "partner_end_bp": PARTNER_END_BP, "partner_short_bp": PARTNER_SHORT_BP,
                                            "depth_min": DEPTH_MIN, "reciprocal_min": RECIPROCAL_MIN,
                                            "support_min_finds": SUPPORT_MIN_FINDS},
                             "interpretation_guard": "Reference-guided completion is homology evidence: a partner is a "
                                                     "candidate missing piece, a split gene is two pieces of one "
                                                     "reference protein at facing contig ends. No contigs are joined, "
                                                     "no product is named, and no score or confidence changes."}
    result["reference_completion"] = block
    pairs = result.get("ranked_pairs") or []
    if tier != "FULL":
        join_pairs(pairs, tier, {}, [], [], {})
        return block
    db = load_mibig_db(mibig_db)
    block["mibig_db_manifest"] = {k: db["manifest"].get(k) for k in ("clusters", "proteins", "files")}
    cores = _core_reference(reference_map, edge_ids, db["genes"])
    source = {k: "knownclusterblast" for k in cores}
    by_id = {b.bgc_id: b for b in bgcs}
    prots, regions = genome_from_zip(zip_path, label, bgcs, edge_ids)
    depth_of = depth_of or (lambda a, b: None)
    tmp = tempfile.TemporaryDirectory(prefix="rggmci_completion_") if work_dir is None else None
    work = Path(tmp.name if tmp else work_dir)
    work.mkdir(parents=True, exist_ok=True)

    def db_search(query_fasta, _db=None, threads=threads, sensitivity=sensitivity, evalue="1e-5", max_target_seqs=25):
        try:
            hs = align(aligner, binary, Path(query_fasta), Path(db["faa"]), work, threads=threads,
                       sensitivity=sensitivity or "default", subject_db=db["dmnd"] if aligner == "diamond" else "",
                       evalue=evalue, max_target_seqs=max_target_seqs)
        except (subprocess.CalledProcessError, OSError) as exc:
            return {"ok": False, "hits": [], "reason": f"{aligner} search failed: {type(exc).__name__}"}
        return {"ok": True, "hits": hs, "reason": f"ok ({aligner})"}

    try:
        # reference discovery for edge regions with no KnownClusterBlast MIBiG hit in the database
        undiscovered = [b for b in bgcs if b.bgc_id in edge_ids and b.bgc_id not in cores]
        block["discovery"] = {"regions": len(undiscovered), "found": 0, "status": "not needed"}
        if undiscovered:
            queries = {b.bgc_id: [k for k, p in prots.items() if p["contig"] == b.contig
                                  and p["start"] < int(b.end or 0) and p["end"] > int(b.start or 0)]
                       for b in undiscovered}
            queries = {k: v for k, v in queries.items() if v}
            try:
                found = discover_references(queries, prots, db["root"], threads=threads, sensitivity=sensitivity,
                                            search=db_search, lengths_kb=db["lengths_kb"] or None,
                                            protein_counts=db["protein_counts"]) if queries else {}
                for bgc_id, d in sorted(found.items()):
                    cores[bgc_id], source[bgc_id] = d["accession"], (
                        f"discovered ({d['proteins']} proteins; size by {d['size_check']}"
                        + (f"; tied with {d['runner_up']}" if d["tied"] else "") + ")")
                block["discovery"].update(found=len(found), status="run")
            except RuntimeError as exc:
                block["discovery"]["status"] = f"not run: {exc}"
        cores = dict(sorted(cores.items()))
        refs_needed = sorted(set(cores.values()))
        qfaa, gfaa = work / "references.faa", work / "genome.faa"
        _write_faa(qfaa, ((g["id"], g["aa"]) for acc in refs_needed for g in db["genes"][acc]))
        _write_faa(gfaa, ((k, v["aa"]) for k, v in prots.items()))
        hits = align(aligner, binary, qfaa, gfaa, work, threads=threads, sensitivity=sensitivity) if refs_needed else []
        by_acc: dict[str, list[dict]] = {}
        for h in hits:
            by_acc.setdefault(mibig_accession(h["qseqid"]), []).append(h)
        staged = []
        for bgc_id, acc in cores.items():
            b = by_id[bgc_id]
            core = {"contig": b.contig, "start": int(b.start or 0), "end": int(b.end or 0), "edge": "True"}
            ref = db["genes"][acc]
            rows, _ = search(ref, prots, regions, by_acc.get(acc, []), core)
            splits, status = split_genes(ref, prots, regions, by_acc.get(acc, []))
            staged.append((bgc_id, acc, core, ref, rows, splits, status))
        # the reciprocal test searches every core's finds against all MIBiG proteins: one search per genome, cached by
        # sequence, instead of one database load per core
        cache: dict[str, list[dict]] = {}
        wanted = sorted({seq for _b, _a, core, ref, rows, _s, _st in staged
                         for seq in _partner_finds(rows, prots, by_acc.get(_a, []), core, ref)[1].values()})

        def cached_search(query_fasta, _db=None, **kw):
            recs, name = [], None
            for line in open(query_fasta):
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    recs.append([name, ""])
                elif recs:
                    recs[-1][1] += line
            missing = sorted({q for _, q in recs if q not in cache})
            if missing:
                mf = work / "reciprocal.faa"
                _write_faa(mf, ((f"m{i}", q) for i, q in enumerate(missing)))
                res = db_search(str(mf), threads=kw.get("threads", threads), sensitivity=kw.get("sensitivity"),
                                max_target_seqs=kw.get("max_target_seqs", 0))
                if not res.get("ok"):
                    return res
                for q in missing:
                    cache[q] = []
                for h in res["hits"]:
                    cache[missing[int(h["qseqid"][1:])]].append(h)
            return {"ok": True, "hits": [dict(h, qseqid=n) for n, q in recs for h in cache.get(q, [])],
                    "reason": f"ok ({aligner}, one search per genome)"}

        if wanted:
            with tempfile.NamedTemporaryFile("w", suffix=".faa", dir=work, delete=False) as fh:
                fh.write("".join(f">w{i}\n{q}\n" for i, q in enumerate(wanted)))
            cached_search(fh.name, sensitivity=sensitivity)
        # the housekeeping check scans the neighbours of lone finds with Pfam: a first pass without Pfam names the lone
        # finds of every core, one Pfam scan covers all their neighbours, and the second pass reads that scan
        pfam_cache: dict[str, list[str]] = {}
        if pfam_hmm:
            lone_pids = set()
            for bgc_id, acc, core, ref, rows, splits, status in staged:
                partner_checks(rows, prots, by_acc.get(acc, []), core, acc, mibig_db=db["root"], threads=threads,
                               sensitivity=sensitivity, search=cached_search, compounds=db["compounds"], ref=ref)
                lone_pids |= {r["best_protein"] for r in rows if r.get("partner_verdict") == "SINGLE_GENE"}
            wanted_n = {k for ks in neighbour_ids(prots, sorted(lone_pids)).values() for k in ks}
            if wanted_n:
                pfam_cache = _pfam_names({k: prots[k]["aa"] for k in sorted(wanted_n)}, str(pfam_hmm), threads,
                                         only=HOUSEKEEPING_PFAM)
                pfam_cache = {k: pfam_cache.get(k, []) for k in wanted_n} if pfam_cache else {}

        def cached_pfam(seqs, hmm, cpus=threads):
            missing = {k: v for k, v in seqs.items() if k not in pfam_cache}
            if missing:
                got = _pfam_names(missing, hmm, cpus, only=HOUSEKEEPING_PFAM)
                if not got:
                    return {}
                pfam_cache.update({k: got.get(k, []) for k in missing})
            return {k: pfam_cache[k] for k in seqs}

        per_core = []
        for bgc_id, acc, core, ref, rows, splits, status in staged:
            checks = partner_checks(rows, prots, by_acc.get(acc, []), core, acc, mibig_db=db["root"],
                                    pfam_hmm=pfam_hmm, threads=threads, sensitivity=sensitivity,
                                    search=cached_search, pfam_names=cached_pfam, compounds=db["compounds"], ref=ref)
            per_core.append((bgc_id, acc, core, ref, rows, splits, status, checks))
        all_splits = []
        for bgc_id, acc, core, ref, rows, splits, status, checks in per_core:
            for s in splits:
                s["reference"], s["_core_bgc"] = acc, bgc_id
                s["_piece_bgcs"] = {(region_of(prots[x], regions) or {}).get("bgc_id") for x in s["_pieces"]} - {None}
            all_splits += splits
        block["recurrent_common_gene_rows"] = mark_recurrent_splits(all_splits, db["compounds"])
        all_partners = []
        for bgc_id, acc, core, ref, rows, splits, status, checks in per_core:
            ident = next(r["identity"] for r in regions if r["bgc_id"] == bgc_id)
            comp = "; ".join(db["compounds"].get(acc, []))
            common = {"completion_tier": tier, "aligner": aligner, "core_identity": ident, "reference": acc,
                      "reference_source": source[bgc_id], "reference_compounds": comp}
            live = [s for s in splits if s["split_call"] != "RECURRENT_COMMON_GENE"]
            parts = partner_tests(core, ref, rows, live, prots, regions, by_acc.get(acc, []), depth_of)
            for p in parts:
                p["_core_bgc"] = bgc_id
            all_partners += parts
            block["cores"].append({**common, "split_check": status, "reference_genes": len(ref),
                                   "present_in_core": sum(r["status"] == "PRESENT_IN_CORE" for r in rows),
                                   "found_clear_elsewhere": sum(r["status"] == "MISSING_FOUND_CLEAR" for r in rows),
                                   "split_genes_clear": sum(s["split_call"] == "CLEAR" for s in splits),
                                   "partner_verdicts": "; ".join(f"{k}:{v}" for k, v in
                                                                 sorted(checks["verdicts"].items())),
                                   "reciprocal_check": checks["reciprocal_mibig"],
                                   "neighbour_pfam_check": checks["neighbour_pfam"],
                                   "accepted_partners": sum(p["accepted"] for p in parts)})
            block["genes"] += [{**common, **r} for r in rows]
            block["splits"] += [{**common, **{k: v for k, v in s.items() if not k.startswith("_")}} for s in splits]
            block["partners"] += [{**common, **{k: v for k, v in p.items() if not k.startswith("_")}} for p in parts]
        join_pairs(pairs, tier, cores, all_splits, all_partners, prots)
        block["genome_proteins"] = len(prots)
        block["references_searched"] = len(refs_needed)
    finally:
        if tmp:
            tmp.cleanup()
    return block


def write_tables(block: dict[str, Any], out_dir: str | Path, prefix: str) -> list[Path]:
    """The three completion tables as TSV (`<prefix>reference_completion.tsv`, `…split_genes.tsv`,
    `…partner_contigs.tsv`). Written even when empty, so a reader sees the tier."""
    out = Path(out_dir)
    paths = []
    for name, cols, key in (("reference_completion", COMPLETION_GENE_COLS, "genes"),
                            ("split_genes", COMPLETION_SPLIT_COLS + ["split_call_before_recurrence"], "splits"),
                            ("partner_contigs", COMPLETION_PARTNER_COLS, "partners")):
        path = out / f"{prefix}{name}.tsv"
        with open(path, "w", newline="") as fh:
            w = SafeWriter(fh, delimiter="\t")
            w.writerow(cols)
            rows = block.get(key) or []
            if not rows:
                w.writerow([block.get("completion_tier", "")] + [""] * (len(cols) - 1))
            for r in rows:
                w.writerow([r.get(c, "") for c in cols])
        paths.append(path)
    return paths


def main(argv: list[str] | None = None) -> int:
    """`python -m mamey.ref_completion build-mibig-db <MIBiG gbk folder> <out folder> [--compounds <index.json>]`."""
    import argparse
    ap = argparse.ArgumentParser(prog="ref_completion", description="Build the MIBiG protein database that RG-GMCI's "
                                 "reference-guided completion reads. Offline: it never downloads anything.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build-mibig-db", help="write a MIBiG protein database from a folder of MIBiG GenBank files")
    b.add_argument("gbk_dir", type=Path, help="folder of MIBiG GenBank files, one cluster per file")
    b.add_argument("out_dir", type=Path, help="a new or empty folder")
    b.add_argument("--compounds", type=Path, help="JSON with an 'entries' list of {accession, compounds}")
    b.add_argument("--diamond", type=Path, help="DIAMOND binary for the .dmnd (default: $RGGMCI_DIAMOND, then PATH)")
    a = ap.parse_args(argv)
    m = build_mibig_db(a.gbk_dir, a.out_dir, compounds_json=a.compounds, diamond=a.diamond)
    sys.stdout.write(f"{m['clusters']} clusters, {m['proteins']} proteins -> {a.out_dir}"
                     + ("" if m["diamond_db"] else " (no DIAMOND found: FASTA only; BLAST+ can still use it)") + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
