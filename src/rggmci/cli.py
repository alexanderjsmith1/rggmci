"""rggmci — surface candidate BGC fragments that may belong to one pathway, from antiSMASH results.

One genome:   rggmci genome.zip --out result.json --pairs pairs.tsv
Many genomes: rggmci zips/ more.zip --out-dir rggmci_results/
BLASTp layer: rggmci fasta genome.zip --out-dir queries/        (then run NCBI BLASTp yourself)
              rggmci blastp-layer --manifest queries/blastp_manifest.json --hits Hit_Table.csv --xml2 r.xml --out-dir layer/
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import run
from .csv_safety import SafeDictWriter

GROUP_COLS = ["genome", "group", "n_regions", "n_contigs", "n_at_contig_ends", "high_pairs", "regions", "pairs",
              "possible_moderate_links"]
SUMMARY_COLS = ["genome", "regions", "scored_pairs", "high", "moderate", "high_cross_contig", "high_both_at_contig_ends",
                "candidate_groups", "largest_group_regions", "related_loci", "reference_map_status",
                "completion_tier", "completion_partners_accepted", "completion_split_genes_clear", "error"]


def _inputs(paths: list[Path]) -> list[Path]:
    out: list[Path] = []
    for p in paths:
        if p.is_dir():
            out += sorted(x for x in p.rglob("*.zip") if not x.name.startswith("._"))
        else:
            out.append(p)
    return out


def _write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = SafeDictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def _summary(name: str, res: dict) -> dict:
    pairs = res.get("ranked_pairs") or []
    conf = [str(p.get("rggmci_confidence", "")) for p in pairs]
    groups = res.get("candidate_groups") or []
    return {"genome": name, "regions": len(res.get("regions", [])), "scored_pairs": len(pairs),
            "high": sum(c.startswith("HIGH") for c in conf), "moderate": sum(c.startswith("MODERATE") for c in conf),
            "high_cross_contig": sum(1 for p in pairs if str(p.get("rggmci_confidence", "")).startswith("HIGH")
                                     and p.get("contig_a") != p.get("contig_b")),
            # A subset of high_cross_contig: two regions at the two ends of ONE contig are not a split pair.
            "high_both_at_contig_ends": sum(1 for p in pairs if str(p.get("rggmci_confidence", "")).startswith("HIGH")
                                            and p.get("contig_a") != p.get("contig_b")
                                            and p.get("both_at_contig_ends")),
            "candidate_groups": len(groups), "largest_group_regions": max((g["n_regions"] for g in groups), default=0),
            "related_loci": len(res.get("related_locus_pairs") or []),
            "reference_map_status": res.get("reference_map_status"),
            "completion_tier": (res.get("reference_completion") or {}).get("completion_tier", ""),
            "completion_partners_accepted": sum(1 for p in (res.get("reference_completion") or {}).get("partners", [])
                                                if p.get("accepted")),
            "completion_split_genes_clear": sum(1 for p in (res.get("reference_completion") or {}).get("splits", [])
                                                if p.get("split_call") == "CLEAR"),
            "error": ""}


def fasta_main(argv) -> int:
    from .blastp import emit_fasta
    ap = argparse.ArgumentParser(prog="rggmci fasta", description="Write BLASTp query FASTA for the regions in "
                                 "RG-GMCI candidate pairs, for you to run on NCBI BLASTp and read back with blastp-layer.")
    ap.add_argument("zip", type=Path, help="one antiSMASH result ZIP")
    ap.add_argument("--out-dir", type=Path, required=True, help="an empty or new folder")
    ap.add_argument("--include-moderate", action="store_true", help="also MODERATE pairs (default: HIGH only)")
    ap.add_argument("--per-region", type=int, default=6, help="biosynthetic proteins per region (default 6)")
    ap.add_argument("--edge-genes", type=int, default=0,
                    help="also this many proteins nearest the contig end, of any kind (default 0)")
    ap.add_argument("--proteins-per-file", type=int, default=20)
    a = ap.parse_args(argv)
    if not a.zip.is_file():
        print(f"not found: {a.zip}", file=sys.stderr)
        return 2
    try:
        m = emit_fasta(a.zip, run(a.zip), a.out_dir, proteins_per_file=a.proteins_per_file,
                       confidences=("HIGH", "MODERATE") if a.include_moderate else ("HIGH",),
                       per_region=a.per_region, edge_genes=a.edge_genes)
    except FileExistsError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(f"{len(m['pairs'])} pairs, {len(m['proteins'])} proteins in {len(m['files'])} FASTA file(s) -> {a.out_dir}. "
          f"See {a.out_dir}/HOW_TO_BLASTP.md")
    return 0


def layer_main(argv) -> int:
    from .blastp import blastp_layer, write_layer
    ap = argparse.ArgumentParser(prog="rggmci blastp-layer", description="Read NCBI BLASTp results for an emitted "
                                 "manifest. The layer sits beside the RG-GMCI confidence and never changes it.")
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--hits", type=Path, nargs="*", default=[], help="NCBI Hit Table (CSV) file(s)")
    ap.add_argument("--xml2", type=Path, nargs="*", default=[], help="NCBI Single-file XML2 file(s): adds organism names")
    ap.add_argument("--exclude-organism", action="append", default=[],
                    help="drop hits whose organism contains this text, e.g. the genome's own species (repeatable)")
    ap.add_argument("--top-n", type=int, default=5, help="organisms kept per protein for the shared-organism check")
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args(argv)
    missing = [str(p) for p in [a.manifest, *a.hits, *a.xml2] if not p.is_file()]
    if missing or not (a.hits or a.xml2):
        print(f"not found: {', '.join(missing)}" if missing else "give --hits and/or --xml2", file=sys.stderr)
        return 2
    layer = blastp_layer(a.manifest, a.hits, a.xml2, top_n=a.top_n, exclude_organisms=a.exclude_organism)
    write_layer(layer, a.out_dir)
    st = [r["blastp_status"] for r in layer["proteins"]]
    print(f"{len(st)} proteins: {st.count('HIT')} with hits, {sum(s.startswith('NO_HITS') for s in st)} no hits, "
          f"{st.count('NOT_IN_RESULTS')} not in the results; {len(layer['pairs'])} pairs -> {a.out_dir}")
    if layer["unmapped_queries"]:
        print(f"{len(layer['unmapped_queries'])} result queries did not match the manifest (renamed headers?); "
              "they were left out", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] in ("--version", "-V"):
        from . import __version__, _provenance
        p = _provenance()
        print(f"rggmci {__version__} (source {p.get('built_from_bundle')}, build {p.get('build')}, "
              f"engine {p.get('engine')})")
        return 0
    if argv and argv[0] == "fasta":
        return fasta_main(argv[1:])
    if argv and argv[0] == "blastp-layer":
        return layer_main(argv[1:])
    if argv and argv[0] == "build-mibig-db":
        from .ref_completion import main as completion_main
        return completion_main(argv)
    ap = argparse.ArgumentParser(prog="rggmci", description="Candidate BGC fragments that may belong to one "
                                 "pathway, from antiSMASH results. Candidates only: check each one at gene level.")
    ap.add_argument("inputs", nargs="+", type=Path, help="antiSMASH result ZIP(s), or folders of them")
    ap.add_argument("--out", type=Path, help="one genome: write the full result as JSON")
    ap.add_argument("--pairs", type=Path, help="one genome: write ranked pairs as TSV")
    ap.add_argument("--out-dir", type=Path, help="many genomes: per-genome JSON, pairs and groups, plus SUMMARY.tsv")
    ap.add_argument("--diamond-db", type=Path, help="optional: a DIAMOND database of MIBiG proteins (ids start with the "
                    "MIBiG accession). Adds residue-tiling evidence to pairs; never changes a confidence")
    ap.add_argument("--diamond", type=Path, help="the DIAMOND binary (default: $RGGMCI_DIAMOND, then diamond on PATH)")
    ap.add_argument("--residue-scope", choices=("st_paralog", "all"), default="st_paralog",
                    help="pairs to test: those demoted by the whole-gene paralog gate (default), or all")
    ap.add_argument("--mibig-db", type=Path, help="MIBiG protein database folder from `rggmci build-mibig-db` "
                    "(default: $RGGMCI_MIBIG_DB). With it and DIAMOND or BLAST+, each edge region's MIBiG reference is "
                    "searched across the whole genome: missing and split genes, and tested partner contigs")
    ap.add_argument("--reference-completion", choices=("auto", "off"), default="auto",
                    help="auto (default) runs reference-guided completion when its inputs are found; off skips it. "
                         "Every pair and table says which tier ran")
    ap.add_argument("--sensitivity", choices=("sensitive", "more-sensitive", "ultra-sensitive", "default"),
                    default="sensitive", help="DIAMOND mode for the completion search (default sensitive)")
    ap.add_argument("--threads", type=int, default=4, help="aligner threads (default 4)")
    ap.add_argument("--pfam", type=Path, help="optional pressed Pfam-A.hmm (needs pyhmmer): sets aside a lone partner "
                    "find beside housekeeping genes")
    a = ap.parse_args(argv)
    residue = {"reference_completion": a.reference_completion, "mibig_db": a.mibig_db, "diamond": a.diamond,
               "sensitivity": a.sensitivity, "threads": a.threads, "pfam_hmm": a.pfam}
    if a.diamond_db:
        residue.update(diamond_db=a.diamond_db, residue_scope=a.residue_scope)
    zips = _inputs(a.inputs)
    missing = [z for z in zips if not z.is_file()]
    if missing or not zips:
        print(f"not found: {', '.join(map(str, missing)) or 'no .zip inputs'}", file=sys.stderr)
        return 2
    if len(zips) > 1 and not a.out_dir:
        print("several inputs: use --out-dir", file=sys.stderr)
        return 2
    if a.out_dir:
        a.out_dir.mkdir(parents=True, exist_ok=True)
    summary, all_groups = [], []
    for z in zips:
        name = z.stem
        try:
            res = run(z, **residue)
        except Exception as exc:   # report and keep going; one bad ZIP should not stop a batch
            summary.append({"genome": name, "error": f"{type(exc).__name__}: {exc}"[:200]})
            print(f"{z.name}: ERROR {type(exc).__name__}", file=sys.stderr)
            continue
        pairs = res.get("ranked_pairs") or []
        s = _summary(name, res)
        summary.append(s)
        for g in res.get("candidate_groups") or []:
            all_groups.append({"genome": name, **{k: g[k] for k in ("group", "n_regions", "n_contigs", "n_at_contig_ends",
                                                                   "high_pairs")},
                               "regions": "; ".join(f"{r['bgc_id']} {r['contig']} [{r['products']}] {r['edge_status']}"
                                                    for r in g["regions"]),
                               "pairs": "; ".join(g["pairs"]),
                               "possible_moderate_links": "; ".join(g["possible_moderate_links"])})
        if a.out:
            a.out.write_text(json.dumps(res, indent=1, default=str) + "\n")
        if a.pairs and pairs:
            _write_tsv(a.pairs, list(pairs[0].keys()), pairs)
        if a.out_dir:
            (a.out_dir / f"{name}.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
            if pairs:
                _write_tsv(a.out_dir / f"{name}_pairs.tsv", list(pairs[0].keys()), pairs)
            related = res.get("related_locus_pairs") or []
            if related:   # interior-region pairs: related loci, never rescues
                _write_tsv(a.out_dir / f"{name}_related_loci.tsv", list(related[0].keys()), related)
            from .ref_completion import write_tables
            write_tables(res.get("reference_completion") or {}, a.out_dir, f"{name}_")
        print(f"{z.name}: {s['regions']} regions, {s['scored_pairs']} pairs, {s['high']} HIGH "
              f"({s['high_cross_contig']} across contigs), {s['candidate_groups']} candidate groups; "
              f"reference completion {s['completion_tier']}")
    if a.out_dir:
        _write_tsv(a.out_dir / "SUMMARY.tsv", SUMMARY_COLS, summary)
        _write_tsv(a.out_dir / "CANDIDATE_GROUPS.tsv", GROUP_COLS, all_groups)
    return 1 if any(s.get("error") for s in summary) else 0
