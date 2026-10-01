"""Lifted by build_rggmci_package.py from Sapote-Mamey source. Do not edit here."""
from __future__ import annotations


# --- from mamey/class_architecture.py: ANTISMASH_PRODUCT_CATEGORY, product_families, product_family
ANTISMASH_PRODUCT_CATEGORY = {
    # nrps
    "cdps": "nrps", "isocyanide-nrp": "nrps", "mycosporine": "nrps", "napaa": "nrps", "nrp-metallophore": "nrps",
    "nrps": "nrps", "nrps-like": "nrps", "t3nrps-iterative": "nrps", "thioamide-nrp": "nrps",
    # pks
    "arylpolyene": "pks", "benzoxazole": "pks", "hgle-ks": "pks", "hr-t2pks": "pks", "pks": "pks",
    "pks-like": "pks", "prodigiosin": "pks", "pufa": "pks", "t1pks": "pks", "t2pks": "pks", "t3pks": "pks",
    "transat-pks": "pks", "transat-pks-like": "pks",
    # ripp
    "atropopeptide": "ripp", "azole-containing-ripp": "ripp", "bottromycin": "ripp", "crocagin": "ripp",
    "cyanobactin": "ripp", "darobactin": "ripp", "fungal-ripp-like": "ripp", "guanidinotides": "ripp",
    "lanthipeptide-class-i": "ripp", "lanthipeptide-class-ii": "ripp", "lanthipeptide-class-iii": "ripp",
    "lanthipeptide-class-iv": "ripp", "lanthipeptide-class-v": "ripp", "lassopeptide": "ripp", "linaridin": "ripp",
    "lipolanthine": "ripp", "methanobactin": "ripp", "proteusin": "ripp", "ranthipeptide": "ripp",
    "redox-cofactor": "ripp", "ripp": "ripp", "ripp-like": "ripp", "rre-containing": "ripp",
    "sactipeptide": "ripp", "thioamitides": "ripp", "triceptide": "ripp",
    # terpene
    "quinone_isoprenoid_chain": "terpene", "terpene": "terpene", "terpene-precursor": "terpene",
    # saccharide
    "oligosaccharide": "saccharide", "saccharide": "saccharide",
    # other
    "2dos": "other", "acyl_amino_acids": "other", "amglyccycl": "other", "aminocoumarin": "other",
    "aminopolycarboxylic-acid": "other", "azoxy-crosslink": "other", "azoxy-dimer": "other",
    "betalactone": "other", "blactam": "other", "butyrolactone": "other", "deazapurine": "other",
    "ectoine": "other", "fatty_acid": "other", "furan": "other", "halogenated": "other", "hserlactone": "other",
    "hydrogen-cyanide": "other", "hydroxytropolone": "other", "indole": "other", "isocyanide": "other",
    "lincosamides": "other", "melanin": "other", "naggn": "other", "ni-siderophore": "other",
    "nucleoside": "other", "opine-like-metallophore": "other", "other": "other", "phenazine": "other",
    "phosphoglycolipid": "other", "phosphonate": "other", "phosphonate-like": "other",
    "polyhalogenated-pyrrole": "other", "polyyne": "other", "pyrrolidine": "other", "resorcinol": "other",
}

def product_family(product) -> str:
    """antiSMASH category of one product type, lower case: nrps, pks, ripp, terpene, saccharide or other.

    A type missing from the table (a newer antiSMASH release) is placed by its name when the name says its
    family ("…-RiPP-like", "…peptide", "…PKS", "…-KS", "…NRP…", "…terpene…"); otherwise it is "other".
    """
    p = str(product).strip().lower().replace("_", "-")
    categories = {k.replace("_", "-"): v for k, v in ANTISMASH_PRODUCT_CATEGORY.items()}
    if p in categories:
        return categories[p]
    if "ripp" in p or "peptide" in p or "lanthi" in p:
        return "ripp"
    if "pks" in p or p.endswith("-ks"):
        return "pks"
    if "nrp" in p:
        return "nrps"
    if "terpene" in p:
        return "terpene"
    if "saccharide" in p:
        return "saccharide"
    return "other"

def product_families(products) -> set[str]:
    """The antiSMASH categories of a region's product types."""
    return {product_family(p) for p in (products or []) if str(p).strip()}
