"""The BGC-record fields the RG-GMCI scorer reads.

The Sapote-Mamey engine's BGCRecord carries many more fields. The scorer in core.py uses only these,
so the standalone package keeps a small record with the same names and the same meanings. The build
checks that every engine BGCRecord field the scorer touches is present here.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class BGCRecord:
    bgc_id: str                      # BGC001 … in the engine's intake order
    contig: str                      # full contig / node id, as antiSMASH wrote it
    region_number: int               # antiSMASH regionNNN
    start: int = 1                   # absolute start on the contig
    end: int = 0                     # absolute end on the contig
    contig_length: int = 0           # true contig length (0 = unknown)
    products: list[str] = field(default_factory=list)
    edge_status: str = "Unknown"     # Interior | Edge | Full-contig | Unknown
    source_gbk: str = ""             # region GBK path inside the ZIP
    subcluster_hits: list = field(default_factory=list)  # filled by the scorer (top-5 SubClusterBlast)
