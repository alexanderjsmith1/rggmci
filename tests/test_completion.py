"""Reference-guided completion in the standalone package: the tier without a database, and the split rules."""
from pathlib import Path

from rggmci.ref_completion import complete, is_mobile_gene, mark_recurrent_splits


def test_without_a_database_every_pair_says_so(tmp_path, monkeypatch):
    monkeypatch.delenv("RGGMCI_MIBIG_DB", raising=False)
    res = {"ranked_pairs": [{"pair": "BGC001+BGC002", "bgc_a": "BGC001", "bgc_b": "BGC002"}]}
    complete(res, Path(tmp_path, "TST-1.zip"), [], {}, set(), mibig_db=None)
    assert res["reference_completion"]["completion_tier"] == "NO_MIBIG_PROTEINS"
    assert res["ranked_pairs"][0]["completion_tier"] == "NO_MIBIG_PROTEINS"


def test_split_rules():
    comp = {"BGC9999001": ["testomycin"], "BGC9999002": ["otheromycin"]}
    real = [{"reference": "BGC9999001", "_core_bgc": "BGC001", "_pieces": ("q1", "q2"),
             "_piece_bgcs": {"BGC001", "BGC002"}, "split_call": "CLEAR"},
            {"reference": "BGC9999002", "_core_bgc": "BGC002", "_pieces": ("q1", "q2"),
             "_piece_bgcs": {"BGC001", "BGC002"}, "split_call": "CLEAR"}]
    assert mark_recurrent_splits(real, comp) == 0
    assert is_mobile_gene({"product": "transposase"}) and not is_mobile_gene({"product": "kinase"})
