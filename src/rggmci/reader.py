"""Read an antiSMASH result ZIP into the records RG-GMCI scores.

Mirrors the part of the Sapote-Mamey engine's `parse_bgcs_from_zip` that sets the six fields the scorer
reads. The helpers it calls (`_parsers.py`, `_ids.py`) are lifted from the engine source at build time,
and a parity test in the engine's own suite compares this reader with the engine on real ZIPs.
"""
from __future__ import annotations

import re
from pathlib import Path

from ._parsers import (_contig_length_map_from_zip, _edge_status, _feature_products, _record_contig_id,
                       _region_orig_bounds_from_zip, _replicon_intake_key, read_genbank_records)
from ._records import BGCRecord


def read_regions(zip_path: str | Path) -> list[BGCRecord]:
    records = read_genbank_records(zip_path, region_only=True) or read_genbank_records(zip_path, region_only=False)
    records.sort(key=_replicon_intake_key)
    orig_bounds = _region_orig_bounds_from_zip(zip_path)
    contig_lengths = _contig_length_map_from_zip(zip_path)
    bgcs: list[BGCRecord] = []
    seen: set = set()
    for name, rec in records:
        is_region_record = "region" in Path(name).name.lower()
        contig_id = _record_contig_id(rec)
        if is_region_record and contig_id not in contig_lengths and rec.id not in contig_lengths:
            contig_len = 0
        else:
            contig_len = contig_lengths.get(contig_id, contig_lengths.get(rec.id, len(rec.seq)))
        is_circular = str((rec.annotations or {}).get("topology", "")).lower() == "circular"
        m = re.search(r"region(\d+)", Path(name).name, flags=re.I)
        region_num = int(m.group(1)) if m else None
        products, start, end = [], 1, contig_len
        if is_region_record and name in orig_bounds:
            start, end = orig_bounds[name]
        for f in rec.features:
            if f.type in {"region", "protocluster", "cand_cluster"}:
                products.extend(_feature_products(f))
                if not (is_region_record and name in orig_bounds):
                    try:
                        start = int(f.location.start) + 1
                        end = int(f.location.end)
                    except (TypeError, ValueError, AttributeError):
                        continue
        products = sorted(set(products))
        if not products:
            for f in rec.features:
                if f.type == "CDS":
                    products.extend(_feature_products(f)[:1])
                    if len(products) >= 3:
                        break
            products = sorted(set(products))
        key = (contig_id, region_num, start, end, tuple(products))
        if key in seen:
            continue
        seen.add(key)
        bgcs.append(BGCRecord(
            bgc_id=f"BGC{len(bgcs) + 1:03d}",
            contig=contig_id,
            region_number=region_num or (len(bgcs) + 1),
            start=start, end=end, contig_length=contig_len,
            products=products,
            edge_status=_edge_status(start, end, contig_len, is_circular=is_circular),
            source_gbk=name,
        ))
    return bgcs
