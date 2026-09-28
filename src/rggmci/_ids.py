"""Lifted by build_rggmci_package.py from Sapote-Mamey source. Do not edit here."""
from __future__ import annotations
import re
import json, re, zipfile, os, statistics, sys
from pathlib import Path


# --- from mamey/crosswalk.py: contig_key
def contig_key(name: str) -> str:
    """Normalise a contig id so the same physical contig matches across files even when the SPAdes coverage
    suffix is formatted differently (e.g. _cov_63.42318 vs _cov_63.042318) — match on node id/length, never
    the cov float (feedback §1 / T-4). Shared util: use for EVERY contig / Source_GBK join, not just orphan
    detection, so cross-file joins (GBK merges, package-dir/strain-id reconciliation) all normalise the same."""
    if not name:
        return ""
    m = re.match(r"(NODE_\d+_length_\d+)", name)
    if m:
        return m.group(1)
    return re.sub(r"[_.]cov[_=].*$", "", name)


# --- from mamey/antismash_evidence.py: _region_key_from_name
def _region_key_from_name(name: str) -> tuple[str | None, int | None, str | None]:
    """Extract a stable antiSMASH region key from a filename.

    antiSMASH restarts the ``_cN``/``regionNNN`` counter on every contig.
    Therefore the bare integer is NOT a valid key for multi-contig assemblies.
    The stable key is the full contig + cluster composite, for example
    ``WOFH01000001.1_c1`` or ``CP073042.1_c38``.
    """
    base = Path(name).name

    # Standard region GBK/TXT names: <contig>.region001.gbk
    m = re.match(r"(.+?)\.region0*(\d+)\.[^.]+$", base, flags=re.I)
    if m:
        contig = m.group(1)
        region = int(m.group(2))
        return contig, region, f"{contig}_c{region}"

    # antiSMASH clusterblast names: <contig>_c1.txt.  Contigs may be accessions
    # like WOFH01000001.1 or CP073042.1, not just NODE_*.
    m = re.match(r"(.+?)_c0*(\d+)\.txt$", base, flags=re.I)
    if m:
        contig = m.group(1)
        region = int(m.group(2))
        return contig, region, f"{contig}_c{region}"

    # Last-resort region number, deliberately not used for KCB assignment unless
    # no contig can be recovered. This prevents cross-contig score pooling.
    m = re.search(r"region0*(\d+)", base, flags=re.I)
    if m:
        return None, int(m.group(1)), None
    return None, None, None
