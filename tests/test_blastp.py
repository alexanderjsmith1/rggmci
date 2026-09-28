"""The BLASTp round trip: emit query FASTA, then read NCBI-shaped results back as a layer. No network."""
import json
import zipfile

import pytest

from rggmci.blastp import blastp_layer, emit_fasta
from rggmci._blastp_io import parse_query_id

AA = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA"


def gbk(name, contig, cds):
    """A minimal antiSMASH-like region GenBank: one region feature and CDS with gene_kind and translation."""
    length = 6000
    feats = [f"     region          1..{length}\n                     /product=\"NRPS\"\n"]
    for i, (kind, start) in enumerate(cds, 1):
        feats.append(f"     CDS             {start}..{start + 599}\n"
                     f"                     /locus_tag=\"{contig}_{i}\"\n"
                     f"                     /gene_kind=\"{kind}\"\n"
                     f"                     /translation=\"{AA}\"\n")
    seq = "a" * length
    origin = "".join(f"{i + 1:>9} {' '.join(seq[i + j:i + j + 10] for j in range(0, 60, 10))}\n"
                     for i in range(0, length, 60))
    return (f"LOCUS       {contig:<16} {length:>11} bp    DNA     linear   BCT 01-JAN-2026\n"
            f"DEFINITION  test.\nACCESSION   {contig}\nVERSION     {contig}\nKEYWORDS    .\nSOURCE      test\n"
            f"  ORGANISM  test\nFEATURES             Location/Qualifiers\n" + "".join(feats) +
            "ORIGIN\n" + origin + "//\n")


@pytest.fixture
def zip_and_result(tmp_path):
    z = tmp_path / "g.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("c1.region001.gbk", gbk("c1.region001", "c1", [("biosynthetic", 100), ("regulatory", 1000),
                                                                   ("biosynthetic-additional", 2000)] +
                                           [("biosynthetic", 2800 + 700 * k) for k in range(4)]))
        zf.writestr("c2.region001.gbk", gbk("c2.region001", "c2", [("biosynthetic", 100), ("other", 900)]))
    result = {"ranked_pairs": [
        {"bgc_a": "BGC001", "bgc_b": "BGC002", "contig_a": "c1", "contig_b": "c2",
         "rggmci_confidence": "HIGH_RG_GMCI_RESCUE"},
        {"bgc_a": "BGC001", "bgc_b": "BGC002", "contig_a": "c1", "contig_b": "c1",
         "rggmci_confidence": "HIGH_RG_GMCI_RESCUE"}]}
    return z, result


def test_emit_takes_biosynthetic_genes_capped_and_headers_parse_back(zip_and_result, tmp_path):
    z, result = zip_and_result
    m = emit_fasta(z, result, tmp_path / "q", per_region=3)
    assert [p["pair"] for p in m["pairs"]] == ["BGC001~BGC002"]          # the same-contig pair is left out
    by = {}
    for p in m["proteins"]:
        by.setdefault(p["bgc_id"], []).append(p)
    assert len(by["BGC001"]) == 3 and all(p["gene_kind"] != "regulatory" for p in by["BGC001"])
    assert [p["gene_kind"] for p in by["BGC002"]] == ["biosynthetic"]
    for p in m["proteins"]:
        meta = parse_query_id(p["query_header"])
        assert meta["bgc_id"] == p["bgc_id"] and meta["gene"] == p["locus_tag"] and meta["node"] == p["node_id"]
    fasta = (tmp_path / "q" / m["files"][0]["file"]).read_text()
    assert fasta.count(">") == len(m["proteins"])
    with pytest.raises(FileExistsError):                                  # earlier queries are never overwritten
        emit_fasta(z, result, tmp_path / "q")


def test_edge_genes_add_the_cds_nearest_the_break(zip_and_result, tmp_path):
    z, result = zip_and_result
    m = emit_fasta(z, result, tmp_path / "q", per_region=1, edge_genes=2)
    roles = [(p["bgc_id"], p["selection_role"]) for p in m["proteins"]]
    assert ("BGC002", "edge") in roles and ("BGC001", "edge") in roles


def _manifest(tmp_path, zip_and_result):
    z, result = zip_and_result
    m = emit_fasta(z, result, tmp_path / "q", per_region=2)
    return tmp_path / "q" / "blastp_manifest.json", m


def _xml(queries):
    """NCBI Single-file XML2 shape: one Search per query; hits as (accession, organism, bitscore)."""
    body = ""
    for qid, qlen, hits in queries:
        hx = "".join(
            f"<Hit><description><HitDescr><id>ref|{acc}|</id><accession>{acc}</accession>"
            f"<title>protein [{org}]</title><sciname>{org}</sciname></HitDescr></description><len>{qlen}</len>"
            f"<hsps><Hsp><bit-score>{bits}</bit-score><evalue>1e-50</evalue><identity>{int(qlen * .8)}</identity>"
            f"<align-len>{qlen}</align-len></Hsp></hsps></Hit>" for acc, org, bits in hits)
        body += (f"<BlastOutput2><report><Report><results><Results><search><Search><query-title>{qid}</query-title>"
                 f"<query-len>{qlen}</query-len><hits>{hx}</hits></Search></search></Results></results></Report>"
                 f"</report></BlastOutput2>")
    return f'<?xml version="1.0"?><BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">{body}</BlastXML2>'


def test_layer_reports_shared_organisms_and_keeps_confidence(tmp_path, zip_and_result):
    mpath, m = _manifest(tmp_path, zip_and_result)
    q = [p["query_header"] for p in m["proteins"]]
    a = [h for h, p in zip(q, m["proteins"]) if p["bgc_id"] == "BGC001"]
    b = [h for h, p in zip(q, m["proteins"]) if p["bgc_id"] == "BGC002"]
    x = tmp_path / "r.xml"
    x.write_text(_xml([(a[0], 190, [("WP_1", "Streptomyces sp. A", 300), ("WP_2", "Nocardia sp. N", 200)]),
                       (a[1], 190, []),
                       (b[0], 190, [("WP_3", "Streptomyces sp. A", 250)])]))
    layer = blastp_layer(mpath, xml2s=[x])
    st = {r["locus_tag"]: r["blastp_status"] for r in layer["proteins"]}
    assert sorted(st.values()) == ["HIT", "HIT", "NO_HITS"]
    pr = layer["pairs"][0]
    assert pr["rggmci_confidence"] == "HIGH_RG_GMCI_RESCUE"                # carried, never changed
    assert pr["organisms_with_homologs_of_both"] == 1 and pr["shared_organisms"] == "Streptomyces sp. A"
    ex = blastp_layer(mpath, xml2s=[x], exclude_organisms=["Streptomyces"])
    assert ex["pairs"][0]["organisms_with_homologs_of_both"] == 0
    assert "NO_HITS_AFTER_EXCLUSION" in {r["blastp_status"] for r in ex["proteins"]}


def test_hit_table_alone_cannot_tell_no_hits_from_not_run(tmp_path, zip_and_result):
    mpath, m = _manifest(tmp_path, zip_and_result)
    first = m["proteins"][0]["query_header"]
    t = tmp_path / "hits.csv"   # headerless, as NCBI writes it
    t.write_text(f"{first},WP_9.1,81.5,190,30,1,1,190,5,194,1e-60,310,90.1\n"
                 f"Query_77,WP_9.1,50,100,50,0,1,100,1,100,1e-5,60,70\n")
    layer = blastp_layer(mpath, hit_tables=[t])
    st = [r["blastp_status"] for r in layer["proteins"]]
    assert st.count("HIT") == 1 and st.count("NOT_IN_RESULTS") == len(st) - 1 and "NO_HITS" not in st
    assert layer["pairs"][0]["shared_organisms"] == "needs XML2"
    assert layer["unmapped_queries"] == ["Query_77"]                        # renamed queries are reported, not guessed
    top = next(r for r in layer["proteins"] if r["blastp_status"] == "HIT")
    assert top["top_identity"] == 81.5
