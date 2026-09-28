"""The package must reproduce the Sapote-Mamey engine's own result on the public example, byte for byte.

`examples/VWPH00000000.1_subset_expected_engine.json` was written by the engine (see its "engine" key) when this repository was generated.
If this test fails, the package and the engine disagree: regenerate the repository from Sapote-Mamey rather
than editing either side here.
"""
import json
from pathlib import Path

import rggmci
from rggmci.reader import read_regions

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ('bgc_id', 'contig', 'region_number', 'contig_length', 'products', 'edge_status')


def test_example_matches_engine():
    z = ROOT / "examples" / "VWPH00000000.1_subset.antismash.zip"
    expected = json.loads((ROOT / "examples" / "VWPH00000000.1_subset_expected_engine.json").read_text())
    b = read_regions(z)
    got = json.loads(json.dumps({"records": [[getattr(x, f) for f in FIELDS] for x in b],
                                 "result": rggmci.run_rggmci(z, b)}, sort_keys=True, default=str))
    assert got["result"]["ranked_pairs"], "the example must score pairs, or this compares nothing"
    assert got["records"] == expected["records"]
    assert got["result"] == expected["result"]
