"""RG-GMCI candidate groups, and one GenBank file per group for the tools that read genes.

RG-GMCI surfaces antiSMASH regions that may belong to one biosynthetic pathway, by biosynthetic logic read through
reference homology. It is not evidence to join contigs: if the reads had supported a join, the assembler would have
made it. Rescues pair only regions antiSMASH flags as on a contig edge, and never two regions on one contig; every
other pair is a related locus and never reaches this module (Alex, 2026-09-27).

1. `candidate_groups`: HIGH pairs on different contigs that share a region, gathered into one group, so a pathway
   spread over several fragments reads as one candidate. MODERATE pairs never join groups: chained together they
   link a genome's similar clusters into one large network. They are listed beside each group as possible further
   members.
2. `write_group_genbanks`: the member regions' own antiSMASH GenBank records, one after another in one file, with a
   members table beside it. Nothing is joined, reordered or filled: each record stays the region antiSMASH wrote.
   The completeness, gene-compare, reference-align and relate tools read every record of a GenBank file, so a group
   file works as one query there.

It never changes a pair's score or confidence. It imports only the standard library and `csv_safety`, which the
standalone `rggmci` package also ships, so the package ships this file unchanged as its `groups.py`.
"""
from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from typing import Any

from .csv_safety import SafeDictWriter

CONTIG_END = {"Edge", "Full-contig"}
GROUP_LABEL = ("RG-GMCI candidate group: fragments on {n} contigs, linked by shared reference homology; "
               "not a contig join.")
MEMBER_COLS = ["order", "bgc_id", "identity", "contig", "antismash_region", "products", "edge_status", "region_file",
               "definition"]


def both_at_contig_ends(pair: dict[str, Any]) -> bool:
    return pair.get("edge_a") in CONTIG_END and pair.get("edge_b") in CONTIG_END


def _cross_contig(pair: dict[str, Any]) -> bool:
    a, b = pair.get("contig_a"), pair.get("contig_b")
    return bool(a) and bool(b) and a != b


def candidate_groups(pairs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Gather HIGH cross-contig pairs that share a region into groups; list touching MODERATE pairs beside them."""
    edges = [p for p in pairs if str(p.get("rggmci_confidence", "")).startswith("HIGH") and _cross_contig(p)]
    moderate = [p for p in pairs if str(p.get("rggmci_confidence", "")).startswith("MODERATE") and _cross_contig(p)]
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for p in edges:
        ra, rb = find(p["bgc_a"]), find(p["bgc_b"])
        if ra != rb:
            parent[ra] = rb
    members: dict[str, dict[str, Any]] = {}
    for p in edges:
        g = members.setdefault(find(p["bgc_a"]), {"regions": {}, "pairs": [], "high_pairs": 0})
        g["regions"][p["bgc_a"]] = (p.get("contig_a"), p.get("products_a"), p.get("edge_a"))
        g["regions"][p["bgc_b"]] = (p.get("contig_b"), p.get("products_b"), p.get("edge_b"))
        g["pairs"].append(p.get("pair") or f"{p['bgc_a']}+{p['bgc_b']}")
        g["high_pairs"] += 1
    out = []
    for g in members.values():
        regs = sorted(g["regions"].items())
        out.append({
            "regions": [{"bgc_id": k, "contig": v[0], "products": v[1], "edge_status": v[2]} for k, v in regs],
            "n_regions": len(regs),
            "n_contigs": len({v[0] for _, v in regs}),
            "n_at_contig_ends": sum(1 for _, v in regs if v[2] in CONTIG_END),
            "pairs": sorted(g["pairs"]),
            "high_pairs": g["high_pairs"],
            "possible_moderate_links": sorted({(p.get("pair") or f"{p['bgc_a']}+{p['bgc_b']}") for p in moderate
                                              if p["bgc_a"] in g["regions"] or p["bgc_b"] in g["regions"]}),
        })
    out.sort(key=lambda g: (-g["high_pairs"], -g["n_regions"]))
    for i, g in enumerate(out, 1):
        g["group"] = f"G{i:02d}"
        g["label"] = GROUP_LABEL.format(n=g["n_contigs"])
    return out


def _definition(text: str) -> str:
    m = re.search(r"^DEFINITION\s+(.*?)(?:\n\S|\Z)", text, re.S | re.M)
    return " ".join(m.group(1).split()) if m else ""


def _member_name(zf: zipfile.ZipFile, source_gbk: str) -> str | None:
    """The archive member for a region file, matched on its path or, failing that, its file name."""
    names = [n for n in zf.namelist() if not n.endswith("/")]
    if source_gbk in names:
        return source_gbk
    base = Path(source_gbk).name
    hits = [n for n in names if Path(n).name == base]
    return hits[0] if len(hits) == 1 else None


def write_group_genbanks(zip_path: str | Path, groups: list[dict[str, Any]], regions: dict[str, dict[str, Any]],
                         out_dir: str | Path, prefix: str) -> list[dict[str, Any]]:
    """Write `<prefix>_<group>.gbk` and `<prefix>_<group>.members.tsv` for each group.

    `regions` maps a region alias to its `strain`, `contig`, `antismash_region`, `source_gbk`, `products` and
    `edge_status`. A group whose region file cannot be found in the archive is reported and not written, so a group
    file always holds every member. Returns one status row per group.
    """
    out_dir = Path(out_dir)
    status = []
    with zipfile.ZipFile(zip_path) as zf:
        for g in groups:
            rows, texts, missing, identity_missing = [], [], [], []
            for i, r in enumerate(g["regions"], 1):
                info = regions.get(r["bgc_id"], {})
                member = _member_name(zf, str(info.get("source_gbk") or ""))
                if member is None:
                    missing.append(r["bgc_id"])
                    continue
                components = [info.get("strain"), info.get("contig"), info.get("antismash_region"), r.get("bgc_id")]
                if not all(isinstance(x, str) and x.strip() and x not in {"None", "?", "IDENTITY_HOLD"}
                           for x in components):
                    identity_missing.append(r["bgc_id"])
                    continue
                text = zf.read(member).decode("utf-8", errors="replace").rstrip("\n")
                texts.append(text if text.endswith("//") else text + "\n//")
                products = info.get("products") or r.get("products") or ""
                rows.append({
                    "order": i, "bgc_id": r["bgc_id"],
                    "identity": " / ".join(str(x) for x in (info.get("strain"), info.get("contig"),
                                                              info.get("antismash_region"), r["bgc_id"])),
                    "contig": info.get("contig", ""), "antismash_region": info.get("antismash_region", ""),
                    "products": "; ".join(products) if isinstance(products, (list, tuple)) else products,
                    "edge_status": info.get("edge_status") or r.get("edge_status") or "",
                    "region_file": member, "definition": _definition(text),
                })
            row = {"group": g["group"], "n_regions": g["n_regions"], "n_contigs": g["n_contigs"], "gbk": "",
                   "members": "", "status": "WRITTEN", "missing_region_files": "; ".join(missing)}
            if missing:
                row["status"] = "NOT_WRITTEN_REGION_FILE_MISSING"
                status.append(row)
                continue
            if identity_missing:
                row["status"] = "NOT_WRITTEN_IDENTITY_HOLD"
                row["identity_holds"] = "; ".join(identity_missing)
                status.append(row)
                continue
            out_dir.mkdir(parents=True, exist_ok=True)
            gbk = out_dir / f"{prefix}_{g['group']}.gbk"
            gbk.write_text("\n".join(texts) + "\n")
            members = out_dir / f"{prefix}_{g['group']}.members.tsv"
            with open(members, "w", newline="") as fh:
                w = SafeDictWriter(fh, fieldnames=MEMBER_COLS, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
            row["gbk"], row["members"] = gbk.name, members.name
            status.append(row)
    return status


def read_members(gbk_path: str | Path) -> list[dict[str, str]]:
    """The members table written beside a group GenBank file, or [] when there is none."""
    p = Path(gbk_path)
    members = p.with_name(p.name[: -len(p.suffix)] + ".members.tsv") if p.suffix else p.with_name(p.name + ".members.tsv")
    if not members.is_file():
        return []
    with open(members, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
