"""A BLASTp round trip for RG-GMCI candidate pairs: emit query FASTA, read the results back as a second layer.

`emit_fasta` writes the proteins of the regions in RG-GMCI candidate pairs as FASTA files a person can paste
into NCBI web BLASTp, plus a manifest. `blastp_layer` reads the Hit Table CSV (and, better, the Single-file
XML2) that NCBI returns, and reports per protein and per pair what came back.

The question the layer helps answer: do close homologs of both fragments' proteins turn up in one source
organism? If they do, that organism may carry an intact version of the pathway, which is worth a look.
It is similarity, not identity, and it never changes a pair's RG-GMCI score or confidence.

The FASTA header, the NCBI-safe batching and the result parsers are lifted from Sapote-Mamey
(`_blastp_io.py`), so results read the same way here and in that project's `blastp-followup`.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
from pathlib import Path
from typing import Any, Iterable

from ._blastp_io import (CLAIM_SAFETY, assign_rounds, best_hits_by_query, fasta_header, merge_hit_xml,
                         parse_hit_table_csv, parse_query_id, parse_xml2, wrap_fasta)
from ._parsers import read_genbank_records
from .csv_safety import SafeDictWriter
from .reader import read_regions

EDGE_FLANK_BP = 5000          # the same flank the reader uses for edge_status
KIND_ORDER = {"biosynthetic": 0, "biosynthetic-additional": 1}
MANIFEST = "blastp_manifest.json"


def _q(feature: Any, key: str) -> str:
    v = (getattr(feature, "qualifiers", None) or {}).get(key)
    return str(v[0]) if v else ""


def _token(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.:+-]+", "_", s).strip("_") or "genome"


def select_proteins(zip_path: str | Path, result: dict[str, Any], *, confidences: Iterable[str] = ("HIGH",),
                    cross_contig_only: bool = True, per_region: int = 6, edge_genes: int = 0,
                    min_aa: int = 30) -> tuple[list[dict], list[dict]]:
    """Return (pairs, protein rows) for the regions in the chosen RG-GMCI pairs.

    Per region: up to `per_region` biosynthetic genes (antiSMASH gene_kind biosynthetic, then
    biosynthetic-additional), taken from the whole region. When the region touches a contig end they are
    ordered nearest-the-end first, so the cap keeps the genes by the break. `edge_genes` adds that many
    CDS nearest the contig end whatever their kind, for checking whether a gene cut by the break continues
    on the other fragment.
    """
    conf = tuple(confidences)
    pairs = [p for p in (result.get("ranked_pairs") or [])
             if str(p.get("rggmci_confidence", "")).startswith(conf)
             and (not cross_contig_only or p.get("contig_a") != p.get("contig_b"))]
    wanted: dict[str, list[str]] = {}
    for p in pairs:
        pid = f'{p["bgc_a"]}~{p["bgc_b"]}'
        wanted.setdefault(p["bgc_a"], []).append(pid)
        wanted.setdefault(p["bgc_b"], []).append(pid)
    if not wanted:
        return pairs, []
    bgcs = {b.bgc_id: b for b in read_regions(zip_path)}
    recs: dict[str, Any] = {}
    for name, rec in (read_genbank_records(zip_path, region_only=True) or read_genbank_records(zip_path)):
        recs.setdefault(name, rec)
    rows: list[dict] = []
    for bid in sorted(wanted):
        b = bgcs[bid]
        rec = recs[b.source_gbk]
        region_file = "region" in Path(b.source_gbk).name.lower()
        offset = b.start - 1 if region_file else 0          # region GBK coordinates start at the region
        lo, hi = (0, len(rec.seq)) if region_file else (b.start - 1, b.end)
        left = b.start <= EDGE_FLANK_BP
        right = bool(b.contig_length) and (b.contig_length - b.end) <= EDGE_FLANK_BP

        def to_end(f: Any) -> int:
            s, e = int(f.location.start) + offset, int(f.location.end) + offset
            d = ([s] if left else []) + ([b.contig_length - e] if right else [])
            return min(d) if d else s

        cds = [f for f in rec.features if f.type == "CDS" and _q(f, "translation")
               and lo <= int(f.location.start) and int(f.location.end) <= hi
               and len(_q(f, "translation")) >= min_aa]
        core = sorted((f for f in cds if _q(f, "gene_kind") in KIND_ORDER),
                      key=lambda f: (KIND_ORDER[_q(f, "gene_kind")], to_end(f)))[:per_region]
        chosen = [(f, "core" if _q(f, "gene_kind") == "biosynthetic" else "additional") for f in core]
        if edge_genes and (left or right):
            have = {id(f) for f, _ in chosen}
            chosen += [(f, "edge") for f in sorted(cds, key=to_end)[:edge_genes] if id(f) not in have]
        for slot, (f, role) in enumerate(chosen, 1):
            seq = _q(f, "translation").rstrip("*")
            rows.append({
                "bgc_id": bid, "slot": slot, "selection_role": role,
                "locus_tag": _q(f, "locus_tag") or _q(f, "protein_id") or _q(f, "gene"),
                "protein_id": _q(f, "protein_id"), "node_id": b.contig,
                "antismash_region": f"region{b.region_number:03d}",
                "start": int(f.location.start) + offset + 1, "end": int(f.location.end) + offset,
                "aa_len": len(seq), "sequence": seq, "gene_kind": _q(f, "gene_kind"),
                "antismash_product": _q(f, "product"), "bp_to_contig_end": to_end(f) if (left or right) else "",
                "selection_reason": f"rggmci_{role}", "pairs": ";".join(wanted[bid]),
            })
    return pairs, rows


def emit_fasta(zip_path: str | Path, result: dict[str, Any], out_dir: str | Path, *, genome: str | None = None,
               proteins_per_file: int = 20, max_residues: int = 85000, **select: Any) -> dict[str, Any]:
    """Write query FASTA files, a manifest (JSON and TSV) and HOW_TO_BLASTP.md into `out_dir`."""
    out = Path(out_dir)
    if (out / MANIFEST).exists():
        raise FileExistsError(f"{out / MANIFEST} exists; choose an empty folder so earlier queries stay intact")
    out.mkdir(parents=True, exist_ok=True)
    genome = _token(genome or Path(zip_path).stem)
    pairs, rows = select_proteins(zip_path, result, **select)
    files = []
    for i, batch in enumerate(assign_rounds(rows, genome, proteins_per_file=proteins_per_file,
                                            max_residues=max_residues), 1):
        name = f"{genome}_blastp_{i:02d}.fasta"
        text = ""
        for r in batch:
            r["query_header"] = fasta_header(genome, r)[1:]
            r["fasta_file"] = name
            text += wrap_fasta(">" + r["query_header"], r["sequence"])
        (out / name).write_text(text)
        files.append({"file": name, "proteins": len(batch), "residues": sum(r["aa_len"] for r in batch)})
    keep = ("rggmci_confidence", "rggmci_score", "contig_a", "contig_b", "edge_a", "edge_b",
            "products_a", "products_b", "both_at_contig_ends")
    manifest = {
        "genome": genome, "antismash_zip": Path(zip_path).name,
        "antismash_zip_sha256": hashlib.sha256(Path(zip_path).read_bytes()).hexdigest(),
        "selection": {k: (list(v) if isinstance(v, tuple) else v) for k, v in select.items()},
        "pairs": [{"pair": f'{p["bgc_a"]}~{p["bgc_b"]}', "bgc_a": p["bgc_a"], "bgc_b": p["bgc_b"],
                   **{k: p.get(k) for k in keep}} for p in pairs],
        "proteins": [{k: v for k, v in r.items() if k != "sequence"} for r in rows],
        "files": files,
    }
    (out / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    cols = ["fasta_file", "query_header", "bgc_id", "node_id", "antismash_region", "locus_tag", "selection_role",
            "gene_kind", "antismash_product", "aa_len", "bp_to_contig_end", "pairs"]
    _tsv(out / "blastp_manifest.tsv", cols, rows)
    (out / "HOW_TO_BLASTP.md").write_text(HOW_TO.format(n=len(rows), f=len(files)))
    return manifest


HOW_TO = """# Running these queries on NCBI BLASTp

{n} proteins in {f} FASTA file(s). Each file stays under NCBI's web limit for one submission.

1. Open https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins and paste one FASTA file (or upload it).
2. Database: nr (non-redundant). ClusteredNR returns one representative per cluster, so organism names
   thin out; use it only when nr is too slow, and say so when you report.
3. If this genome is already in NCBI, exclude its own organism in "Organism → exclude", or its own proteins
   come back as 100% hits. You can also pass `--exclude-organism` when reading the results.
4. When the search finishes, open "Download All" and save both:
   - **Hit Table (CSV)**: every hit, with identity and scores;
   - **Single-file XML2**: the same hits with organism names, and it lists queries that found nothing.
     Without it, a protein missing from the Hit Table could mean "no hits" or "not run".
5. Read them back:

       rggmci blastp-layer --manifest blastp_manifest.json --hits Hit_Table.csv --xml2 results.xml --out-dir layer/

Do not rename the query headers; they carry the region, gene and pair each result belongs to.
BLASTp results are sequence similarity. They do not show that two fragments are one pathway.
"""


def _tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = SafeDictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})


def _organism(hit: dict[str, Any]) -> str:
    name = str(hit.get("subject_sciname") or "").strip()
    if not name:   # XML2 without a sciname: take the last [..] of the title, as deposited
        m = re.findall(r"\[([^\[\]]+)\]", str(hit.get("subject_title") or ""))
        name = m[-1].strip() if m else ""
    return name


def blastp_layer(manifest_path: str | Path, hit_tables: Iterable[str | Path] = (), xml2s: Iterable[str | Path] = (),
                 *, top_n: int = 5, exclude_organisms: Iterable[str] = ()) -> dict[str, Any]:
    """Read NCBI results for an emitted manifest. Returns protein rows, pair rows and unmapped queries."""
    m = json.loads(Path(manifest_path).read_text())
    excl = [e.lower() for e in exclude_organisms if e.strip()]
    kept = lambda org: bool(org) and not any(e in org.lower() for e in excl)
    hits = [h for t in hit_tables for h in parse_hit_table_csv(t)]
    xml: dict[str, dict[str, Any]] = {}
    for x in xml2s:
        xml.update(parse_xml2(x))
    if hits and xml:
        hits = merge_hit_xml(hits, xml)
    key = lambda qid: (parse_query_id(qid).get("bgc_id"), parse_query_id(qid).get("gene"))
    wanted = {(p["bgc_id"], parse_query_id(p["query_header"]).get("gene")) for p in m["proteins"]}
    table_best = {key(q): r for q, r in best_hits_by_query(hits).items()}
    xml_by = {key(q): v for q, v in xml.items()}
    unmapped = sorted({q for q in [h.query_id for h in hits] + list(xml) if key(q) not in wanted})

    prot_rows: list[dict] = []
    for p in m["proteins"]:
        k = (p["bgc_id"], parse_query_id(p["query_header"]).get("gene"))
        xq = xml_by.get(k)
        xhits = sorted(xq["hits"], key=lambda h: -(h.get("bitscore") or 0)) if xq else []
        if excl:
            xhits = [h for h in xhits if kept(_organism(h))]
        top_orgs = list(dict.fromkeys(o for o in (_organism(h) for h in xhits) if o))[:top_n]
        best = table_best.get(k)
        if best is not None and excl and xq is not None:
            best = None   # the table row may be an excluded organism; take the top kept XML2 hit instead
        if xq is not None:
            status = "HIT" if xhits else ("NO_HITS_AFTER_EXCLUSION" if xq["hits"] else "NO_HITS")
        else:
            status = "HIT" if k in table_best else "NOT_IN_RESULTS"
        row = {**{c: p.get(c) for c in ("bgc_id", "node_id", "antismash_region", "locus_tag", "selection_role",
                                        "gene_kind", "antismash_product", "aa_len", "pairs")},
               "blastp_status": status, "top_organisms": "; ".join(top_orgs)}
        if best is not None:
            row.update(top_identity=round(best.pct_identity, 1), top_query_coverage=(
                "" if best.query_coverage is None else round(best.query_coverage, 3)),
                top_evalue=best.evalue, top_subject=best.subject_accession or best.subject_id,
                top_title=best.subject_title, top_organism=best.subject_sciname)
        elif xhits:
            h = xhits[0]
            qlen = xq.get("query_len") or p.get("aa_len")
            row.update(top_identity=(round(100 * h["identity"] / h["align_len"], 1)
                                     if h.get("identity") and h.get("align_len") else ""),
                       top_query_coverage=round(min(1.0, h["align_len"] / qlen), 3) if h.get("align_len") and qlen else "",
                       top_evalue=h.get("evalue", ""), top_subject=h.get("subject_accession") or h.get("subject_id"),
                       top_title=h.get("subject_title", ""), top_organism=_organism(h))
        row["_orgs"] = set(top_orgs)
        prot_rows.append(row)

    by_bgc: dict[str, list[dict]] = {}
    for r in prot_rows:
        by_bgc.setdefault(r["bgc_id"], []).append(r)
    pair_rows = []
    for pr in m["pairs"]:
        out = {"pair": pr["pair"], "rggmci_confidence": pr.get("rggmci_confidence"),
               "contig_a": pr.get("contig_a"), "contig_b": pr.get("contig_b")}
        orgs = {}
        for side, bid in (("a", pr["bgc_a"]), ("b", pr["bgc_b"])):
            rs = by_bgc.get(bid, [])
            ids = [r["top_identity"] for r in rs if isinstance(r.get("top_identity"), (int, float))]
            out[f"{side}_proteins"] = len(rs)
            out[f"{side}_with_hits"] = sum(r["blastp_status"] == "HIT" for r in rs)
            out[f"{side}_no_hits"] = sum(r["blastp_status"].startswith("NO_HITS") for r in rs)
            out[f"{side}_not_in_results"] = sum(r["blastp_status"] == "NOT_IN_RESULTS" for r in rs)
            out[f"{side}_median_top_identity"] = round(statistics.median(ids), 1) if ids else ""
            orgs[side] = set().union(*[r["_orgs"] for r in rs]) if rs else set()
        have_xml = any(r["top_organisms"] for r in by_bgc.get(pr["bgc_a"], []) + by_bgc.get(pr["bgc_b"], []))
        shared = sorted(orgs["a"] & orgs["b"])
        out["organisms_with_homologs_of_both"] = len(shared) if have_xml else ""
        out["shared_organisms"] = "; ".join(shared[:10]) if have_xml else "needs XML2"
        pair_rows.append(out)
    for r in prot_rows:
        r.pop("_orgs")
    return {"proteins": prot_rows, "pairs": pair_rows, "unmapped_queries": unmapped,
            "claim_safety": CLAIM_SAFETY, "excluded_organisms": list(exclude_organisms)}


def write_layer(layer: dict[str, Any], out_dir: str | Path) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pc = ["bgc_id", "node_id", "antismash_region", "locus_tag", "selection_role", "gene_kind", "antismash_product",
          "aa_len", "blastp_status", "top_identity", "top_query_coverage", "top_evalue", "top_subject", "top_organism",
          "top_title", "top_organisms", "pairs"]
    _tsv(out / "blastp_proteins.tsv", pc, layer["proteins"])
    qc = ["pair", "rggmci_confidence", "contig_a", "contig_b"] + [f"{s}_{c}" for s in "ab" for c in (
        "proteins", "with_hits", "no_hits", "not_in_results", "median_top_identity")] + [
        "organisms_with_homologs_of_both", "shared_organisms"]
    _tsv(out / "blastp_pairs.tsv", qc, layer["pairs"])
    (out / "blastp_layer.json").write_text(json.dumps(layer, indent=1, default=str) + "\n")
