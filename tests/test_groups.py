from rggmci.groups import both_at_contig_ends, candidate_groups


def pair(a, b, ca, cb, ea, eb, conf="HIGH_RG_GMCI_RESCUE"):
    return {"bgc_a": a, "bgc_b": b, "contig_a": ca, "contig_b": cb, "edge_a": ea, "edge_b": eb,
            "products_a": "NRPS", "products_b": "NRPS", "rggmci_confidence": conf, "pair": f"{a}+{b}"}


def test_contig_end_description():
    assert both_at_contig_ends(pair("B1", "B2", "c1", "c2", "Edge", "Full-contig"))
    assert not both_at_contig_ends(pair("B1", "B2", "c1", "c2", "Edge", "Interior"))


def test_interior_high_pairs_are_described_not_hidden():
    g = candidate_groups([pair("B1", "B2", "c1", "c2", "Interior", "Edge")])
    assert len(g) == 1 and g[0]["n_regions"] == 2 and g[0]["n_at_contig_ends"] == 1


def test_same_contig_pairs_do_not_form_groups():
    assert candidate_groups([pair("B1", "B2", "c1", "c1", "Edge", "Edge")]) == []


def test_moderate_pairs_are_listed_beside_groups_not_joined():
    ps = [pair("B1", "B2", "c1", "c2", "Edge", "Full-contig"),
          pair("B2", "B3", "c2", "c3", "Full-contig", "Edge", conf="MODERATE_RG_GMCI_CANDIDATE"),
          pair("B4", "B5", "c4", "c5", "Edge", "Edge", conf="LOW_SHARED_REFERENCE_SIGNAL")]
    g = candidate_groups(ps)
    assert len(g) == 1
    assert [r["bgc_id"] for r in g[0]["regions"]] == ["B1", "B2"]
    assert g[0]["possible_moderate_links"] == ["B2+B3"]


def test_moderate_chains_do_not_build_a_network():
    ps = [pair(f"B{i}", f"B{i+1}", f"c{i}", f"c{i+1}", "Edge", "Edge", conf="MODERATE_RG_GMCI_CANDIDATE")
          for i in range(10)]
    assert candidate_groups(ps) == []


def test_result_names_its_input_and_build(tmp_path):
    import zipfile
    from rggmci import _input_record
    z = tmp_path / "g.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("x.txt", "x")
    rec = _input_record(z)
    assert rec["antismash_zip"] == "g.zip" and len(rec["antismash_zip_sha256"]) == 64
    assert "rggmci_built_from_bundle" in rec


def test_version_is_single_sourced():
    import rggmci
    from rggmci.cli import main
    assert rggmci.__version__ not in ("", "unknown")
    assert main(["--version"]) == 0
