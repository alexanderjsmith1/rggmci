"""The depth fields and the residue-tiling placement, on synthetic ids."""
from rggmci.core import contig_depth, depth_fields
from rggmci.residue_tiling import compare_sides, place_side


def test_depth():
    assert contig_depth("NODE_1_length_9_cov_40.5") == 40.5
    assert contig_depth("TST-1_ctg1") is None
    assert depth_fields("NODE_1_length_9_cov_90", "NODE_2_length_9_cov_30")["depth_flag"] == "DEPTH_MISMATCH"


def test_split_gene_versus_copy():
    a = place_side([("TST-A|1", "BGC9999999|1", 60.0, 1, 1200, 1, 1200, 1500.0)])
    b = place_side([("TST-B|1", "BGC9999999|1", 60.0, 1, 1300, 1300, 2600, 1600.0)])
    c = place_side([("TST-C|1", "BGC9999999|1", 60.0, 1, 1200, 1, 1200, 1400.0)])
    assert compare_sides(a, b)["class"] == "COMPLEMENTARY"
    assert compare_sides(a, c)["class"] == "OVERLAPPING"
