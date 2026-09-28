"""Reference-Guided Genome Mining Candidate Inference (RG-GMCI).

Full, pre-triage linkage pass for fragmented antiSMASH assemblies.

RG-GMCI asks whether two or more fragmented BGC regions in the query genome
map to overlapping, adjacent, or nearby segments of the same reference producer
cluster/genome in antiSMASH ClusterBlast / KnownClusterBlast output.  It is a
homology-guided candidate-linkage method, not a physical contig joiner.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any, Iterable
import itertools
import math
from ._ids import contig_key
import re
import zipfile

from ._records import BGCRecord
from ._ids import _region_key_from_name  # lifted verbatim from mamey/antismash_evidence.py
from .ziputil import regular_file_names
from ._families import product_families  # lifted verbatim from mamey/class_architecture.py


RGGMCI_SCORE_HIGH = 14
RGGMCI_SCORE_MODERATE = 9
OVERLAP_FRACTION_MIN = 0.20  # min overlap/segment-span to count as true OVERLAPPING (else ADJACENT)
RGGMCI_MAX_PER_BGC_REFS = 40
RGGMCI_MAX_EVIDENCE_ROWS = 5000
RGGMCI_MAX_RANKED_PAIRS = 1000
# Locus-number proxy (v9.7.39, P-4'): when antiSMASH emits no subject coordinates (the subject
# gene-location table is empty in this format), fall back to subject locus-tag NUMBERS as a coarse
# position proxy so adjacency is computable. Tagged adjacency_basis=LOCUS_PROXY — never a coordinate claim.
# v9.7.41: ADJ_MAX_LOCUS_GAP tightened 60->30 from the 5-strain PRIVATE calibration
# (5 internal PRIVATE strains, ~394 ADJACENT proxy pairs): genuine adjacency is gap~0 (median 0 on all
# five), gap<=30 retains ~95% of ADJACENT pairs, and the 31-60 band is a thin tail (~11 pairs total).
# ADJ_MAX_SPAN kept at 400: max real ADJACENT span across the five strains was 315-400.
ADJ_MAX_LOCUS_GAP = 30
ADJ_MAX_SPAN = 400

# v9.7.42: component-degree guard. With Bug A/B fixed, MIBiG cluster references count toward geometry, and
# on a heavily fragmented assembly a promiscuous hub fragment can become pairwise-HIGH with many partners
# (recalibration: one VERY_POOR strain showed 183 HIGH across 33 BGCs at degree 8-19 — a dense web, not 183 real co-clusters). A
# genuine cross-contig split rarely produces more than ~4-5 confidently-linked fragments (a 2-fragment
# core+tailoring split is degree 1), so a BGC co-clustering with MORE than this is treated as a hub and its
# HIGH pairs demote to MODERATE. Calibrated on the 10-strain set: the clean degree-1 indolocarbazole
# survives; the dense hub webs on the most fragmented strains collapse.
RGGMCI_MAX_HUB_DEGREE = 4

# ── Subject-gene tiling (v9.7.100, P-ST) ─────────────────────────────────────────────────
# The ranked-pairs summary scored pairs on adjacency geometry + identity + product class, but never
# exposed *which subject genes of the shared reference each contig hits*. That subject-gene map is the
# evidence that actually distinguishes a genuine cross-contig cluster split (the two query contigs hit
# DISJOINT, COMPLEMENTARY runs of the same reference cluster — core fragment vs tailoring fragment) from
# shared-machinery paralogy (the two contigs hit the SAME subject genes — same conserved enzymes, two
# independent loci). It was computable from data already parsed onto ClusterBlastReference.subjects, but
# was dropped before output, so a split call rested on the score alone and needed manual KCB-TXT reading
# to confirm. This block computes the tiling per shared reference and aggregates it per pair. It is
# OBSERVATIONAL ONLY in this cut — it adds evidence fields, it does not change score or confidence.
# A minimum of this many subject hits on EACH side is required before a reference is allowed to report a
# COMPLEMENTARY_DISJOINT tiling verdict (one-gene-vs-one-gene is too thin to call complementary).
ST_MIN_SUBJECTS_PER_SIDE = 2
# Pair-level corroboration thresholds: a split verdict requires disjoint complementary tiling on at least
# this many independent shared references, and this many aggregated distinct subject genes per side, so a
# single thin disjoint reference (the default for two unrelated BGCs) is NOT mistaken for a split.
ST_MIN_DISJOINT_REFS = 2
ST_MIN_AGG_SUBJECTS_PER_SIDE = 3
# A split verdict is only surfaced where adjacency geometry already supports linkage on at least this many
# shared references (good_geometry = OVERLAPPING/ADJACENT segments); otherwise disjoint subjects are just
# the default for unrelated BGCs and the pair is reported INSUFFICIENT, not a split.
ST_MIN_GEOMETRY_FOR_SPLIT = 1
# A reference is still complementary when the two fragments share a few of its genes but cover different genes
# otherwise: one gene family present in two copies within one cluster (sugar genes, RRE/methyltransferases) is hit
# from both sides. Two copies of the same machinery share far more (a desferrioxamine pair shares 25-45%).
ST_MAX_SHARED_SUBJECTS = 2
ST_MAX_SHARED_FRACTION = 0.15

# ── Terminus-truncation rescue (v9.7.100, P-TT) ──────────────────────────────────────────
# The simplest, most certain rescue is physical, not homological: an Edge region whose boundary IS the
# contig terminus is a cluster sliced by the assembly break, and its missing arm is on another contig. RG-
# GMCI was going straight to ClusterBlast homology + subject tiling and never checking the coordinate that
# proves truncation. Worse, when the severed arm happened to share a gene class that legitimately occurs in
# multiple copies within ONE cluster (e.g. the Bottromycin_Methyltransferase_RRE — bottromycin clusters
# encode several RRE/methyltransferase genes), subject tiling saw a "shared subject" and called the pair
# OVERLAPPING_PARALOG, vetoing a real split. This block (a) detects terminus truncation from start/end vs
# contig_length, and (b) lets a coordinate-confirmed truncation + same-class partner OVERRIDE a paralogy
# verdict that rests only on a multi-copy-within-cluster gene. The paralog-overlap function is NOT removed —
# it still fires for genuine paralogs (distinct contigs, both core-complete, no terminus truncation).
TT_TERMINUS_MARGIN_BP = 500   # an Edge region whose start<=margin or end>=contig_len-margin is terminus-truncated
TT_SMALL_PARTNER_BP = 25000   # a Full-contig partner this short is a candidate broken-off arm
TT_COMPLEXITY_WINDOW_BP = 200
TT_HOMOPOLYMER_RUN_BP = 20
TT_MIN_SHANNON_BITS = 1.20


def _terminus_truncation(bgc: BGCRecord) -> str:
    """Return 'start', 'end', 'both', or '' for how a region abuts its contig terminus.

    Only meaningful for Edge regions on linear contigs; a region flush with a contig boundary is truncated.
    """
    if getattr(bgc, "edge_status", "") != "Edge":
        return ""
    s = getattr(bgc, "start", None); e = getattr(bgc, "end", None)
    L = getattr(bgc, "contig_length", None)
    if s is None or e is None or not L:
        return ""
    at_start = s <= TT_TERMINUS_MARGIN_BP
    at_end = (L - e) <= TT_TERMINUS_MARGIN_BP
    if at_start and at_end:
        return "both"
    return "start" if at_start else ("end" if at_end else "")


def _terminus_sequence_state(bgc: BGCRecord, contigs: dict[str, str] | None) -> str:
    """Classify sequence complexity at the exact truncated terminus.

    Low-complexity or homopolymer sequence cannot establish a unique assembly
    junction. Such a terminus is routed to ``LONG_READ_ONLY`` downstream; this
    detector never guesses the missing arm.
    """
    side = _terminus_truncation(bgc)
    if not side or not contigs:
        return "NOT_APPLICABLE" if not side else "TERMINUS_SEQUENCE_UNAVAILABLE"
    seq = contigs.get(bgc.contig)
    if seq is None:
        matches = [value for key, value in contigs.items() if contig_key(key) == contig_key(bgc.contig)]
        if len(matches) != 1:
            return "TERMINUS_SEQUENCE_UNAVAILABLE"
        seq = matches[0]
    seq = re.sub(r"[^ACGT]", "", str(seq).upper())
    if not seq:
        return "TERMINUS_SEQUENCE_UNAVAILABLE"
    windows = []
    if side in ("start", "both"):
        windows.append(seq[:TT_COMPLEXITY_WINDOW_BP])
    if side in ("end", "both"):
        windows.append(seq[-TT_COMPLEXITY_WINDOW_BP:])
    for window in windows:
        if re.search(rf"([ACGT])\1{{{TT_HOMOPOLYMER_RUN_BP - 1},}}", window):
            return "HOMOPOLYMER_TERMINUS"
        total = len(window)
        entropy = -sum(
            (window.count(base) / total) * math.log2(window.count(base) / total)
            for base in "ACGT" if window.count(base)
        )
        if entropy < TT_MIN_SHANNON_BITS:
            return "LOW_COMPLEXITY_TERMINUS"
    return "COMPLEX_TERMINUS"


_CLASS_TOKENS = (
    "bottromycin", "lanthipeptide", "lassopeptide", "lap", "thiopeptide", "sactipeptide",
    "linaridin", "lanthidin", "indolocarbazole", "enediyne", "siderophore",
    "ni-siderophore", "nrp-metallophore", "metallophore", "azole-containing-ripp",
    "bacteriocin", "lipopeptide", "glycopeptide", "aminoglycoside",
)
# Generic umbrella tokens that are too common to anchor a terminus rescue on their own (the existing P-2
# gate treats these the same way). A terminus rescue needs either a SPECIFIC shared class from _CLASS_TOKENS
# or the small-severed-arm partner geometry — never a bare saccharide/pks/nrps/ripp match.
_TT_GENERIC_TOKENS = {"saccharide", "pks", "t1pks", "t2pks", "t3pks", "nrps", "ripp",
                      "rre-containing", "terpene", "other", "fatty_acid"}


def _label_class_tokens(products: str) -> set[str]:
    """Map each discrete product label to its single LONGEST-matching _CLASS_TOKENS entry.

    v9.7.377 (AUDIT correctness fix): naive substring containment over the whole
    products string let a QUALIFIED class collide with its own unqualified substring — e.g.
    'NI-siderophore' contains the substring 'siderophore', and 'NRP-metallophore' contains
    'metallophore'. Both qualified/unqualified pairs are separately listed in _CLASS_TOKENS
    because antiSMASH treats them as DIFFERENT, biosynthetically distinct product types (an
    NRPS-independent siderophore/metallophore pathway is not the same machinery as the
    generic NRPS-dependent one) — yet the old substring check credited both sides with the
    unqualified token whenever only the qualified one was present, so
    _shared_class_tokens('NI-siderophore', 'siderophore') falsely returned {'siderophore'}.
    That false "same specific class" verdict is exactly the biosynthetic-logic-complementarity
    proof this rescue is gated on, so the collision could flip a genuine OVERLAPPING_PARALOG
    verdict (a real shared subject gene, e.g. a conserved siderophore-uptake receptor) into a
    false TERMINUS_TRUNCATION_SPLIT claim. Splitting each product string into discrete labels
    first, then keeping only the LONGEST _CLASS_TOKENS substring per label, means a qualified
    label is credited with its qualified token ONLY — 'lanthipeptide-class-i' still correctly
    maps to 'lanthipeptide' (no other token collides there), but 'ni-siderophore' now maps to
    'ni-siderophore' only, never also to 'siderophore'.
    """
    labels = [t.strip().lower() for t in re.split(r"[;,/|]+", products or "") if t.strip()]
    out: set[str] = set()
    for label in labels:
        hits = [tok for tok in _CLASS_TOKENS if tok in label]
        if hits:
            out.add(max(hits, key=len))
    return out


def _shared_class_tokens(products_a: str, products_b: str) -> set[str]:
    """Specific shared biosynthetic class tokens (excludes generic umbrellas)."""
    return _label_class_tokens(products_a) & _label_class_tokens(products_b)


def _subject_tiling(a_subjects: Iterable[str], b_subjects: Iterable[str]) -> dict[str, Any]:
    """Compare the subject-gene sets two query contigs hit in ONE shared reference cluster.

    Returns counts and a per-reference tiling_class:
      COMPLEMENTARY_DISJOINT  — both sides hit >= ST_MIN_SUBJECTS_PER_SIDE subjects and the sets are
                                disjoint, or share at most ST_MAX_SHARED_SUBJECTS genes that are at most
                                ST_MAX_SHARED_FRACTION of all genes hit, each side still hitting
                                ST_MIN_SUBJECTS_PER_SIDE of its own (the split fingerprint: each contig tiles a
                                different part of the same reference cluster).
      OVERLAPPING_SUBJECTS    — the sides share more subject genes than that (paralogy / shared machinery;
                                NOT a clean split — the contigs compete for the same reference genes).
      SINGLETON_OR_THIN       — at least one side hits fewer than ST_MIN_SUBJECTS_PER_SIDE subjects, so
                                complementarity cannot be called either way.
    Subject identity is compared on the bare locus tag as parsed (e.g. 'ABC02802.1', 'QJU69504.1').
    """
    a_set = {s for s in a_subjects if s}
    b_set = {s for s in b_subjects if s}
    shared = a_set & b_set
    a_only = a_set - b_set
    b_only = b_set - a_set
    few_shared = (len(shared) <= ST_MAX_SHARED_SUBJECTS
                  and len(shared) <= ST_MAX_SHARED_FRACTION * len(a_set | b_set)
                  and len(a_only) >= ST_MIN_SUBJECTS_PER_SIDE and len(b_only) >= ST_MIN_SUBJECTS_PER_SIDE)
    if len(a_set) < ST_MIN_SUBJECTS_PER_SIDE or len(b_set) < ST_MIN_SUBJECTS_PER_SIDE:
        tiling = "SINGLETON_OR_THIN"
    elif shared and not few_shared:
        tiling = "OVERLAPPING_SUBJECTS"
    else:
        tiling = "COMPLEMENTARY_DISJOINT"
    return {
        "subjects_a": tuple(sorted(a_set)),
        "subjects_b": tuple(sorted(b_set)),
        "shared_subjects": tuple(sorted(shared)),
        "a_only_subjects": tuple(sorted(a_only)),
        "b_only_subjects": tuple(sorted(b_only)),
        "n_shared_subjects": len(shared),
        "n_a_only_subjects": len(a_only),
        "n_b_only_subjects": len(b_only),
        "subject_tiling_class": tiling,
    }

# v9.7.41 P-2: product-class compatibility gate on ACCEPT. The 5-strain calibration showed gg=1
# single-reference over-promotion as the dominant precision leak; gg>=2 (in _geometry_gate) catches
# most of it, but a residue of pairs reach HIGH on gg>=2 whose only shared "product" is a
# permanent-exclusion token (NAPAA/saccharide) or "other" -- e.g. a both-Interior pair (BGC014+BGC040,
# shared token = NAPAA only). This gate requires a specific shared product class OR a known compatible
# hybrid for HIGH, AFTER excluding other/saccharide/NAPAA and the generic RiPP umbrella. Pairs with no
# specific class on either side (product-dark) are NOT promoted to HIGH on product grounds -> demoted
# HIGH->MODERATE. Demote-only; never promotes. (A permissive-on-dark variant -- keep product-dark pairs
# HIGH-eligible on geometry alone -- is a one-line change if dark-BGC HIGH-rescue is wanted; see field note.)
_P2_EXCL_TOKENS = {"other", "saccharide", "aminopolycarboxylic-acid", "napaa"}
_P2_RIPP = {"lassopeptide", "lanthipeptide", "lanthipeptide-class-i", "lanthipeptide-class-ii",
            "lanthipeptide-class-iii", "lanthipeptide-class-iv", "lanthipeptide-class-v",
            "thiopeptide", "ripp-like", "sactipeptide", "lap", "ras-ripp", "lassopeptide-class",
            "azole-containing-ripp", "ranthipeptide", "thioamitide", "lanthidin", "epipeptide",
            "cyclic-lactone-autoinducer", "spliceotide", "crocagin", "darobactin"}
_P2_UMBRELLA = {"ripp"}
_P2_HYBRID_OK = {frozenset({"nrps", "pks"}), frozenset({"nrps", "t1pks"}),
                 frozenset({"nrps-like", "t1pks"}), frozenset({"pks", "nrps-like"}),
                 frozenset({"nrps", "nrps-like"}), frozenset({"t1pks", "t2pks"})}


def _p2_specific_classes(products: str) -> set[str]:
    """Lower-cased specific product classes, minus permanent exclusions and the RiPP umbrella."""
    cls = {p.strip().lower() for p in re.split(r"[;,]", products or "") if p.strip()}
    # Family words count too, so an NRPS fragment and a PKS-family fragment (arylpolyene, PKS-like, T3PKS) form the
    # NRPS+PKS hybrid below, as they did when products carried antiSMASH's category words.
    cls |= product_families(cls)
    return {c for c in cls if c not in _P2_EXCL_TOKENS and c not in _P2_UMBRELLA}



# ── R1/R3 helper constants (v9.7.61) ────────────────────────────────────────────────────
# R1: tokens that mark a BGC as a standing-rule DROP class (saccharide/NAPAA/hglE-KS).
# When BOTH sides of a HIGH pair are exclusively DROP-class tokens, cap at MODERATE.
# Note: 'pks' alone does not rescue a pair — hglE-KS glycolipid loci have pks as their
# only backbone, so pks+hgle-ks combined is still a DROP-class pair.
_R1_DROP_ONLY_TOKENS = {"other", "saccharide", "napaa", "hgle-ks", "hgle", "pks", "pks-like"}

# R3: tokens that are informationally vacuous for co-cluster assessment.
# When BOTH sides contain only these tokens (no specialist class on either side),
# the pair lacks the class signal needed to confirm a genuine co-cluster → cap at MODERATE.
# Rationale: fatty_acid×fatty_acid and other×other pairs are dominated by conserved
# housekeeping/lipid biosynthesis convergence, not split pathway rescue.
_R3_NOISE_TOKENS = {"other", "saccharide", "napaa", "hgle-ks", "fatty_acid"}


def _is_drop_only(products: str) -> bool:
    """True iff all product tokens in `products` are R1 DROP-class tokens."""
    cls = {t.strip().lower() for t in re.split(r"[;,/|]+", products or "") if t.strip()}
    return bool(cls) and cls <= _R1_DROP_ONLY_TOKENS


def _is_noise_only(products: str) -> bool:
    """True iff all product tokens in `products` are R3 noise tokens (no specialist class)."""
    cls = {t.strip().lower() for t in re.split(r"[;,/|]+", products or "") if t.strip()}
    return bool(cls) and cls <= _R3_NOISE_TOKENS

def _product_class_gate(confidence: str, products_a: str, products_b: str) -> tuple[str, str]:
    """v9.7.41 P-2. Demote HIGH->MODERATE unless the pair has a specific shared product class or a
    known compatible hybrid (after excluding other/saccharide/NAPAA + the RiPP umbrella). Demote-only."""
    if confidence != "HIGH_RG_GMCI_RESCUE":
        return confidence, "OK"
    sa = _p2_specific_classes(products_a)
    sb = _p2_specific_classes(products_b)
    if not sa or not sb:
        return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_no_specific_product_class_after_exclusions"
    # two DIFFERENT specific RiPP subclasses are distinct terminal machineries -> not one cluster
    if sa <= _P2_RIPP and sb <= _P2_RIPP and not (sa & sb):
        return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_distinct_ripp_subclasses"
    if sa & sb:
        return confidence, "OK_shared_specific_product_class"
    for x in sa:
        for y in sb:
            if frozenset({x, y}) in _P2_HYBRID_OK:
                return confidence, "OK_compatible_hybrid"
    return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_incompatible_product_classes"


RGGMCI_MAX_COCLUSTER_MIBIG_REFS = 6  # v9.7.49 promiscuity cap: above this, shared MIBiG refs signal a
# conserved cassette matching many characterized clusters, not a specific co-cluster -> exemption withheld.


def _cocluster_exempt(product_gate_reason: str, mibig_good_geometry_references: int) -> bool:
    """v9.7.42 — co-cluster carve-out for the product gate, from the v9.7.41 lateral review.

    A genuine cross-contig split is *complementary by construction*: the core fragment and the tailoring
    fragment carry DIFFERENT product labels (indole core + halogenated/saccharide arm; NRPS + PKS). The
    v9.7.41 product gate, built to kill convergence noise, demotes exactly that complementarity as
    `incompatible_product_classes` — so it fights the signal RG-GMCI exists to find. This carve-out
    restores such a pair to its post-geometry grade ONLY when it is corroborated by >=2 OVERLAPPING/ADJACENT
    **MIBiG cluster** references (`BGC#######`) — i.e. both fragments tile the same characterized cluster.

    Critically it fires ONLY on the `incompatible_product_classes` reason (two real, different specific
    classes = complementarity). It never fires on `no_specific_product_class_after_exclusions`
    (convergence on an excluded/empty token — the saccharide/NAPAA case) or on `distinct_ripp_subclasses`
    (mutually exclusive terminal machineries). So saccharide-on-saccharide / NAPAA-on-NAPAA convergence and
    distinct-RiPP pairs stay demoted, and the permanent exclusions are preserved.

    v9.7.49 PROMISCUITY CAP: the exemption is also WITHHELD when the shared-MIBiG count exceeds
    RGGMCI_MAX_COCLUSTER_MIBIG_REFS. A genuine split tiles ONE characterized cluster (a few references); a
    conserved cassette (e.g. an aromatic type-II PKS) tiles MANY distinct MIBiG clusters (~dozens), which is
    promiscuity, not co-location. This is orthogonal to the hub-degree guard — that guard meters fragment-pair
    degree and misses the MIBiG-co-cluster path, where pairwise degree stays low (hub_deg None) even as the
    shared-MIBiG count is high. The two guards together close both over-promotion paths."""
    return ("incompatible_product_classes" in (product_gate_reason or "")
            and 2 <= mibig_good_geometry_references <= RGGMCI_MAX_COCLUSTER_MIBIG_REFS)




def _r1_drop_pair_gate(confidence: str, products_a: str, products_b: str) -> tuple[str, str]:
    """R1 (v9.7.61): If BOTH BGCs are exclusively DROP-class tokens (saccharide/NAPAA/hglE-KS
    backbone), cap HIGH→MODERATE. Pairs of genuinely distinct DROP-class BGCs that share
    convergent lipid/glycan references should not occupy HIGH slots; they belong in MODERATE
    for review. Demote-only; does not affect MODERATE or LOW."""
    if confidence != "HIGH_RG_GMCI_RESCUE":
        return confidence, "OK"
    if _is_drop_only(products_a) and _is_drop_only(products_b):
        return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_R1_both_drop_class"
    return confidence, "OK"


def _r2_strong_ref_gate(confidence: str, strong_supporting_references: int) -> tuple[str, str]:
    """R2 (v9.7.61): Require strong_supporting_references ≥ 1 as a co-condition for HIGH,
    alongside the existing good_geometry_references ≥ 2 requirement.

    strong_supporting_references counts references where BOTH sides of the pair independently
    scored HIGH_RG_GMCI_RESCUE on the per-reference scoring (i.e. per-reference confidence was
    HIGH, not just MODERATE or LOW). A pair with good geometry but zero strong references means
    the geometry came from references that only weakly supported either side individually —
    consistent with conserved domain convergence rather than genuine co-cluster rescue.

    Pairs demoted by R2 remain in MODERATE tier and are not suppressed."""
    if confidence != "HIGH_RG_GMCI_RESCUE":
        return confidence, "OK"
    if strong_supporting_references < 1:
        return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_R2_no_strong_references"
    return confidence, "OK"


def _r3_noise_class_gate(confidence: str, products_a: str, products_b: str) -> tuple[str, str]:
    """R3 (v9.7.61): If NEITHER side carries any specialist class token, cap HIGH→MODERATE.
    Pairs of fatty_acid×fatty_acid, other×other, etc. (exclusively noise-class tokens on both
    sides) lack the biosynthetic class signal needed to confirm a genuine co-cluster rescue —
    even with good geometry and multiple strong references, these pairs reflect conserved
    housekeeping or lipid biosynthesis convergence, not secondary metabolite co-localization.

    Note: hglE-KS pairs are caught first by R1; R3 provides redundant cover plus catches
    fatty_acid and other combinations not covered by R1."""
    if confidence != "HIGH_RG_GMCI_RESCUE":
        return confidence, "OK"
    if _is_noise_only(products_a) and _is_noise_only(products_b):
        return "MODERATE_RG_GMCI_CANDIDATE", "DEMOTED_HIGH_TO_MODERATE_R3_noise_class_both_sides"
    return confidence, "OK"

# Sugar-arm restore (v9.7.444). A fragment antiSMASH labels only "saccharide" can be the deoxysugar tailoring arm of a
# glycosylated product whose backbone sits on another contig (a glycosylated angucycline: the core on one contig,
# the aminosugar cassette on another). The product-class gate sets "saccharide" aside, so such a pair could never hold
# HIGH. The backbone side must carry a glycosyltransferase gene of its own: the enzyme that attaches the sugar.
_SUGAR_ONLY_TOKENS = {"saccharide", "oligosaccharide", "other"}
_SUGAR_ARM_BACKBONE_FAMILIES = {"nrps", "pks", "ripp", "terpene"}
_SUGAR_ARM_OTHER_BACKBONES = {"indole", "nucleoside", "betalactone", "phenazine", "aminocoumarin", "2dos", "amglyccycl"}
_GLYCOSYLTRANSFERASE_RE = re.compile(r"glycosyl ?transferase|glycosyltransf|Glycos_transf|Glyco_trans|UDPGT|SMCOG1062",
                                     re.I)


def _types(products: Iterable[str]) -> set[str]:
    return {str(p).strip().lower() for p in products if str(p).strip()}


def _is_sugar_only(products: Iterable[str]) -> bool:
    t = _types(products)
    return bool(t) and t <= _SUGAR_ONLY_TOKENS and bool(t & {"saccharide", "oligosaccharide"})


def _is_sugar_arm_backbone(products: Iterable[str]) -> bool:
    t = _types(products)
    return bool(product_families(t) & _SUGAR_ARM_BACKBONE_FAMILIES) or bool(t & _SUGAR_ARM_OTHER_BACKBONES)


def _sugar_arm_pair(bgc_a: BGCRecord, bgc_b: BGCRecord, glycosyltransferase_bgcs: set[str] | None) -> bool:
    """One side sugar-only, the other a backbone class whose own region carries a glycosyltransferase gene."""
    for arm, core in ((bgc_a, bgc_b), (bgc_b, bgc_a)):
        if (_is_sugar_only(arm.products) and _is_sugar_arm_backbone(core.products)
                and core.bgc_id in (glycosyltransferase_bgcs or ())):
            return True
    return False


def _glycosyltransferase_bgcs(zip_path: str | Path, bgcs: list[BGCRecord]) -> set[str]:
    """BGC ids whose own region GenBank file names a glycosyltransferase gene (annotation text only)."""
    want = {b.source_gbk: b.bgc_id for b in bgcs if getattr(b, "source_gbk", "")}
    out: set[str] = set()
    if not want:
        return out
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = set(regular_file_names(zf))
            for name, bgc_id in want.items():
                if name in names and _GLYCOSYLTRANSFERASE_RE.search(_read_text(zf, name)):
                    out.add(bgc_id)
    except (OSError, zipfile.BadZipFile):
        return set()
    return out


_EDGE_STATUSES = {"Edge", "Full-contig"}
_CONTIG_EDGE_RE = re.compile(r'/contig_edge="(True|False)"')


def _edge_status_fallback(bgcs: list[BGCRecord]) -> set[str]:
    """Regions the engine places at a contig end (within 5 kb of one, or covering the contig)."""
    return {b.bgc_id for b in bgcs if getattr(b, "edge_status", "") in _EDGE_STATUSES}


def _contig_edge_bgcs(zip_path: str | Path, bgcs: list[BGCRecord]) -> set[str]:
    """Regions antiSMASH itself flags as touching a contig edge (`/contig_edge="True"` on the region feature).

    RG-GMCI pairs only these (Alex, 2026-09-27): a region antiSMASH places inside a contig is whole on that contig,
    and an assembly break cannot have separated it from anything. A region file without the flag falls back to the
    engine's own edge status.
    """
    fallback = _edge_status_fallback(bgcs)
    out: set[str] = set()
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = set(regular_file_names(zf))
            for b in bgcs:
                name = getattr(b, "source_gbk", "")
                m = _CONTIG_EDGE_RE.search(_read_text(zf, name)) if name in names else None
                if m:
                    if m.group(1) == "True":
                        out.add(b.bgc_id)
                elif b.bgc_id in fallback:
                    out.add(b.bgc_id)
    except (OSError, zipfile.BadZipFile):
        return fallback
    return out


def _apply_hub_degree_guard(ranked: list[dict[str, Any]]) -> None:
    """v9.7.42 component-degree guard (graph-level, runs last; mutates `ranked` in place).

    With MIBiG references counting toward geometry, a promiscuous hub fragment on a fragmented assembly can
    become pairwise-HIGH with many partners, inflating the lead board (one strain: 183 HIGH across 33 BGCs at
    degree 8-19 — a dense web, not 183 real co-clusters). A genuine cross-contig split rarely produces more
    than RGGMCI_MAX_HUB_DEGREE confidently-linked fragments (a 2-fragment split is degree 1), so a BGC whose
    HIGH degree exceeds the cap is a hub: its HIGH pairs demote HIGH->MODERATE. Each pair records
    `max_endpoint_hub_degree`. Degree is computed once on the post-gate HIGH set (single pass, deterministic);
    a clean low-degree co-cluster (e.g. a degree-1 indolocarbazole) survives, a hub web collapses to MODERATE."""
    degree: dict[str, int] = {}
    for r in ranked:
        if r["rggmci_confidence"] == "HIGH_RG_GMCI_RESCUE":
            degree[r["bgc_a"]] = degree.get(r["bgc_a"], 0) + 1
            degree[r["bgc_b"]] = degree.get(r["bgc_b"], 0) + 1
    for r in ranked:
        d = max(degree.get(r["bgc_a"], 0), degree.get(r["bgc_b"], 0))
        r["max_endpoint_hub_degree"] = d
        if r["rggmci_confidence"] == "HIGH_RG_GMCI_RESCUE" and d > RGGMCI_MAX_HUB_DEGREE:
            r["rggmci_confidence"] = "MODERATE_RG_GMCI_CANDIDATE"
            reason = f"DEMOTED_HUB_PROMISCUITY_degree_{d}_gt_{RGGMCI_MAX_HUB_DEGREE}"
            r["acceptance_gate"] = reason if r["acceptance_gate"] == "OK" else f"{r['acceptance_gate']}+{reason}"


ACCESSION_RE = re.compile(
    # v9.7.42 (Bug A): every branch now requires a digit, so all-caps class/header keywords
    # (NRPS, BLAST, PKS, ...) no longer match as "accessions"; MIBiG BGC####### is matched explicitly
    # (it was previously unmatched, routing MIBiG cluster references to a garbage key).
    r"\b(?:"
    r"BGC\d{6,7}"                                              # MIBiG cluster accession
    r"|(?:NZ_|NC_|NG_|NT_|NW_|NM_|NR_)[A-Z]{0,2}\d[A-Z0-9_]*"  # RefSeq
    r"|C[PM]\d{4,}"                                            # GenBank complete genome (CP/CM)
    r"|[A-Z]{4}\d{2,}"                                         # 4-letter WGS master
    r"|[A-Z]{1,2}\d{5,6}"                                      # v9.7.166: classic GenBank (e.g. DQ149987) — subclusterblast refs use this older format
    r")(?:\.\d+)?\b"
)


@dataclass(frozen=True)
class ClusterBlastReference:
    bgc_id: str
    contig: str
    region_number: int | None
    region_key: str | None
    ref: str
    source: str
    reference_type: str
    rank: int
    nprot: int | None
    cumulative_score: float | None
    mean_identity: float | None
    interval_start: int | None
    interval_end: int | None
    source_file: str
    subjects: tuple[str, ...] = ()   # subject locus tags from the Blast-hits table (P-4' proxy input)
    db_kind: str = "clusterblast"    # P-CBDB v9.7.100: which DB this block came from —
                                     # "knownclusterblast" (characterized MIBiG), "clusterblast" (cross-genome
                                     # GenBank neighbors), or "subclusterblast" (sub-operon; excluded from
                                     # cluster geometry). Previously all three were blind-merged because the
                                     # substring "clusterblast" matches all of them.

    @property
    def interval(self) -> tuple[int | None, int | None]:
        return (self.interval_start, self.interval_end)


def _read_text(zf: zipfile.ZipFile, name: str) -> str:
    return zf.read(name).decode("utf-8", errors="replace")


def _norm_type_tokens(text: str) -> set[str]:
    text = (text or "").lower().replace("_", "-")
    toks = re.split(r"[^a-z0-9+-]+", text)
    return {t for t in toks if len(t) >= 3 and t not in {"type", "other", "unknown", "cluster", "biosynthetic"}}


def _product_tokens(products: Iterable[str]) -> set[str]:
    # antiSMASH's family words count beside the types ("pks" for T3PKS or arylpolyene), as the category words did when
    # products carried them: two fragments of one PKS pathway can be typed differently.
    products = list(products)
    return _norm_type_tokens(";".join(products + sorted(product_families(products))))


def _first_accession(text: str) -> str | None:
    m = ACCESSION_RE.search(text or "")
    return m.group(0) if m else None


def _extract_source(block: str) -> str:
    # antiSMASH has used several header labels across versions.
    for pat in [r"Source\s*:\s*(.+)", r"Organism\s*:\s*(.+)", r"Description\s*:\s*(.+)"]:
        m = re.search(pat, block, flags=re.I)
        if m:
            return re.sub(r"\s+", " ", m.group(1).strip())[:240]
    first = next((ln.strip() for ln in block.splitlines() if ln.strip()), "")
    return re.sub(r"\s+", " ", first)[:240]


def _extract_type(block: str) -> str:
    for pat in [r"Type\s*:\s*(.+)", r"BGC type\s*:\s*(.+)", r"Product\s*:\s*(.+)"]:
        m = re.search(pat, block, flags=re.I)
        if m:
            return re.sub(r"\s+", " ", m.group(1).strip())[:200]
    return ""


def _extract_score(block: str) -> float | None:
    for pat in [r"Cumulative\s+BLAST\s+score\s*:\s*([0-9][0-9,]*(?:\.\d+)?)", r"\bscore\s*[:=]\s*([0-9][0-9,]*(?:\.\d+)?)"]:
        m = re.search(pat, block, flags=re.I)
        if m:
            try:
                return float(m.group(1).replace(",", ""))
            except ValueError:
                return None
    return None


def _extract_nprot(block: str) -> int | None:
    pats = [
        r"Number\s+of\s+proteins\s+with\s+BLAST\s+hits\s+to\s+this\s+cluster\s*:\s*(\d+)",
        r"proteins\s+with\s+BLAST\s+hits[^:\n]*:\s*(\d+)",
        r"\bnprot\s*[:=]\s*(\d+)",
    ]
    for pat in pats:
        m = re.search(pat, block, flags=re.I)
        if m:
            return int(m.group(1))
    return None


def _blast_hits(block: str) -> list[tuple[str, float | None]]:
    """Parse the 'Table of Blast hits' rows: (query, subject, %identity, score, coverage, evalue).

    Returns (subject_locus, %identity) per row. The %identity is column 3 as a BARE number in this
    antiSMASH format (no '%' sign), which is why the old identity-by-'%'-sign scan returned nothing.
    """
    out: list[tuple[str, float | None]] = []
    m = re.search(r"Table of Blast hits.*?\n(.*)", block, flags=re.S)
    if not m:
        return out
    for line in m.group(1).splitlines():
        parts = line.split("\t")
        # v9.7.43: a Blast-hits data row has 6 tab columns (query, subject, %identity, score, %coverage,
        # e-value) with a bare numeric %identity in col 3. The old filter required the query gene to start
        # with 'ctg', which silently dropped every row on genomes that name query genes by locus tag
        # (e.g. GTY48_01375 on public WGS assemblies) -> 0 subjects/identity -> RG-GMCI degenerate. Gate on
        # the row structure instead of the query-gene prefix.
        if len(parts) < 6:
            continue
        try:
            ident = float(parts[2])
        except ValueError:
            continue
        subj = parts[1].strip()
        if subj:
            out.append((subj, ident))
    return out


def _extract_mean_identity(block: str) -> float | None:
    # Primary: mean of the Blast-hits %identity column (col 3, bare numbers in this format).
    idents = [i for _, i in _blast_hits(block) if i is not None]
    if idents:
        return round(sum(idents) / len(idents), 1)
    # Fallback (older formats that DO carry a '%' sign).
    vals: list[float] = []
    for m in re.finditer(r"(?:identity|identities|id)\D{0,20}([0-9]+(?:\.[0-9]+)?)\s*%", block, flags=re.I):
        vals.append(float(m.group(1)))
    if not vals:
        for m in re.finditer(r"([0-9]+(?:\.[0-9]+)?)\s*%", block):
            v = float(m.group(1))
            if 20 <= v <= 100:
                vals.append(v)
    if not vals:
        return None
    return round(sum(vals) / len(vals), 1)


def _extract_interval(block: str) -> tuple[int | None, int | None]:
    """Best-effort reference interval extraction.

    antiSMASH TXT formats differ by version.  We gather plausible coordinate
    spans from reference/subject-side labels and fall back to any coordinate-like
    spans in the block.  The result is only used for coarse overlap/adjacency
    classes, not nucleotide-level claims.
    """
    spans: list[tuple[int, int]] = []
    # v9.7.40 (Fix-1): the Blast-hits rows (query gene starts with 'ctg') carry %identity/%coverage
    # columns whose decimal values (e.g. %coverage 101.5625) were being misparsed as a coordinate span
    # (101-5625), creating phantom intervals that drove false HIGH promotions. clusterblast TXT carries
    # NO reference coordinates; exclude those rows, and never treat a single '.' as a coordinate separator
    # (real ranges use '-' or GenBank '..'). Adjacency on this format comes from the locus-tag proxy, not here.
    scan = "\n".join(l for l in block.splitlines() if not l.strip().startswith("ctg"))
    SEP = r"\s*(?:[-\u2013\u2014]|\.\.)\s*"  # hyphen / en / em dash / GenBank '..' — NOT a lone decimal point

    labeled_patterns = [
        r"(?:Reference|Subject|Hit|Cluster)\s*(?:location|coordinates|range|interval)?[^\n]{0,80}?([0-9]{2,9})" + SEP + r"([0-9]{2,9})",
        r"(?:from|start)\s+([0-9]{2,9})\s+(?:to|end)\s+([0-9]{2,9})",
    ]
    for pat in labeled_patterns:
        for m in re.finditer(pat, scan, flags=re.I):
            a, b = int(m.group(1)), int(m.group(2))
            if a != b:
                spans.append((min(a, b), max(a, b)))

    # GenBank-style accession:start-end strings are especially useful.
    for m in re.finditer(r"(?:\b[A-Z]{2,4}_[A-Z0-9]+(?:\.\d+)?|\bNZ_[A-Z0-9]+(?:\.\d+)?)[:_]([0-9]{2,9})" + SEP + r"([0-9]{2,9})", scan):
        a, b = int(m.group(1)), int(m.group(2))
        if a != b:
            spans.append((min(a, b), max(a, b)))

    # Fall back only if nothing labelled was found.  Filter out tiny spans and
    # very large scaffold-wide coordinates that are unlikely to be locus spans.
    if not spans:
        for m in re.finditer(r"\b([0-9]{3,9})" + SEP + r"([0-9]{3,9})\b", scan):
            a, b = int(m.group(1)), int(m.group(2))
            lo, hi = min(a, b), max(a, b)
            if 20 <= hi - lo <= 250_000:
                spans.append((lo, hi))

    if not spans:
        return None, None
    return min(s[0] for s in spans), max(s[1] for s in spans)


def _split_subject_blocks(text: str) -> list[str]:
    # KnownClusterBlast blocks usually contain one Cumulative BLAST score each.
    parts = re.split(r"\n\s*>>\s*\n|\n\s*>\s*(?=\S)|\n\s*Cluster\s+\d+\s*[:.]", text, flags=re.I)
    blocks = [p for p in parts if re.search(r"Cumulative\s+BLAST\s+score|Number\s+of\s+proteins\s+with\s+BLAST\s+hits|\bType\s*:", p, flags=re.I)]
    return blocks or [text]


def _classify_db_kind(name: str) -> str:
    """Classify a ClusterBlast TXT by its source directory.

    Order matters: 'knownclusterblast' and 'subclusterblast' both CONTAIN the substring 'clusterblast', so
    they must be tested first. Returns 'knownclusterblast', 'subclusterblast', or 'clusterblast'.
    """
    low = name.lower()
    if "knownclusterblast" in low:
        return "knownclusterblast"
    if "subclusterblast" in low:
        return "subclusterblast"
    return "clusterblast"


def _parse_clusterblast_refs_from_text(name: str, text: str, bgc_by_key: dict[str, BGCRecord]) -> list[ClusterBlastReference]:
    contig, region, region_key = _region_key_from_name(name)
    bgc = bgc_by_key.get(region_key or "")
    if not bgc:
        # Region without a stable BGC mapping.  Do not emit loose RG-GMCI refs;
        # pairwise linkage needs locked BGC IDs.
        return []

    refs: list[ClusterBlastReference] = []
    seen: set[tuple[str, int, int | None, int | None]] = set()
    db_kind = _classify_db_kind(name)
    for idx, block in enumerate(_split_subject_blocks(text), 1):
        ref = _first_accession(block)
        if not ref:
            continue
        start, end = _extract_interval(block)
        nprot = _extract_nprot(block)
        score = _extract_score(block)
        ident = _extract_mean_identity(block)
        subjects = tuple(s for s, _ in _blast_hits(block))
        src = _extract_source(block)
        typ = _extract_type(block)
        key = (ref, idx, start, end)
        if key in seen:
            continue
        seen.add(key)
        refs.append(ClusterBlastReference(
            bgc_id=bgc.bgc_id,
            contig=bgc.contig,
            region_number=bgc.region_number,
            region_key=region_key,
            ref=ref,
            source=src,
            reference_type=typ,
            rank=idx,
            nprot=nprot,
            cumulative_score=score,
            mean_identity=ident,
            interval_start=start,
            interval_end=end,
            source_file=name,
            subjects=subjects,
            db_kind=db_kind,
        ))
        if len(refs) >= RGGMCI_MAX_PER_BGC_REFS:
            break
    return refs


def parse_clusterblast_reference_map(zip_path: str | Path, bgcs: list[BGCRecord]) -> dict[str, Any]:
    bgc_by_key = {f"{b.contig}_c{int(b.region_number)}": b for b in bgcs if b.contig and b.region_number is not None}
    refs: list[ClusterBlastReference] = []
    files: list[str] = []
    with zipfile.ZipFile(zip_path) as zf:
        for name in regular_file_names(zf):
            low = name.lower()
            if not low.endswith(".txt"):
                continue
            if "clusterblast" not in low and "knownclusterblast" not in low:
                continue
            files.append(name)
            try:
                refs.extend(_parse_clusterblast_refs_from_text(name, _read_text(zf, name), bgc_by_key))
            except Exception as exc:  # pragma: no cover - defensive parse ledger
                refs.append(ClusterBlastReference(
                    bgc_id="PARSE_ERROR", contig="", region_number=None, region_key=None,
                    ref="PARSE_ERROR", source=f"{name}: {type(exc).__name__}: {exc}", reference_type="",
                    rank=9999, nprot=None, cumulative_score=None, mean_identity=None,
                    interval_start=None, interval_end=None, source_file=name,
                ))
    usable = [r for r in refs if r.bgc_id != "PARSE_ERROR"]
    return {
        "status": "PASS" if usable else "NULL_NO_CLUSTERBLAST_REFERENCES_PARSED",
        "clusterblast_txt_files": files,
        "reference_records": [asdict(r) for r in usable],
        "reference_record_count": len(usable),
        "parse_error_count": len(refs) - len(usable),
        "claim_safety": "ClusterBlast TXT-derived reference map; coordinates are coarse text-parse evidence only.",
    }


def _locus_num(name: str) -> int | None:
    # v9.7.42 (Bug B): strip a trailing GenBank version suffix (.1) so a protein accession like
    # 'QJU69504.1' yields 69504, not the version number 1.
    s = re.sub(r"\.\d+$", "", name or "")
    m = re.findall(r"(\d+)", s)
    return int(m[-1]) if m else None


def _locus_key(name: str):
    """Split a locus tag into (prefix, number): 'AMK09_RS30055' -> ('AMK09_RS', 30055).
    v9.7.42 (Bug B): a GenBank version suffix is stripped first, so 'QJU69504.1' -> ('QJU', 69504)
    rather than ('QJU69504.', 1) — without this, MIBiG/GenBank-protein subjects landed in disjoint
    prefix namespaces and could never form proxy adjacency."""
    s = re.sub(r"\.\d+$", "", name or "")
    m = re.match(r"(.*?)(\d+)$", s)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def _proxy_adjacency(a: ClusterBlastReference, b: ClusterBlastReference):
    """Locus-number proxy used only when antiSMASH emitted no subject coordinates.

    Compares locus NUMBERS only WITHIN a shared locus-tag prefix (namespace) — cross-namespace
    comparison (e.g. AMK09_RS30055 vs B082_RS0106330) is meaningless and would inflate the span.
    Returns (adjacency_class, gap, span), or (None, None, None) when no shared-prefix numeric loci.
    """
    from collections import defaultdict
    pa: dict[str, set[int]] = defaultdict(set)
    pb: dict[str, set[int]] = defaultdict(set)
    for s in (a.subjects or ()):
        pre, num = _locus_key(s)
        if num is not None:
            pa[pre].add(num)
    for s in (b.subjects or ()):
        pre, num = _locus_key(s)
        if num is not None:
            pb[pre].add(num)
    shared = set(pa) & set(pb)
    if not shared:
        return None, None, None
    # The namespace with the smallest gap wins. On a tie, the one with more loci (more evidence) wins, then the
    # prefix. Iterating the set alone made the tie depend on Python's per-process hash seed (v9.7.444).
    ranked = []
    for pre in shared:
        na, nb = pa[pre], pb[pre]
        ranked.append((min(abs(x - y) for x in na for y in nb), -len(na | nb), pre, max(na | nb) - min(na | nb)))
    gap, _, _, span = min(ranked)
    if gap <= ADJ_MAX_LOCUS_GAP and span <= ADJ_MAX_SPAN:
        return "ADJACENT_OR_NEARBY_REFERENCE_SEGMENTS", gap, span
    return "DISTANT_ON_REFERENCE_CAUTION", gap, span


def _adjacency(a: ClusterBlastReference, b: ClusterBlastReference):
    """Returns (adjacency_class, gap, overlap, overlap_fraction, adjacency_basis, span).

    adjacency_basis ∈ {COORDINATE, LOCUS_PROXY, NONE}. LOCUS_PROXY is a coarse position estimate from
    locus-tag numbers (this antiSMASH format emits no subject coordinates) — never a nucleotide-level claim.
    For LOCUS_PROXY, gap/span are locus-number distances; for COORDINATE, gap/overlap are bp.
    """
    a0, a1 = a.interval
    b0, b1 = b.interval
    if None in (a0, a1, b0, b1):
        cls, gap, span = _proxy_adjacency(a, b)
        if cls is not None:
            return cls, gap, None, None, "LOCUS_PROXY", span
        return "SHARED_REFERENCE_NO_INTERVAL", None, None, None, "NONE", None
    assert a0 is not None and a1 is not None and b0 is not None and b1 is not None
    overlap = max(0, min(a1, b1) - max(a0, b0))
    span = min(a1 - a0, b1 - b0)
    ofrac = round(overlap / span, 3) if (span and span > 0) else 0.0
    # true OVERLAPPING only when the overlap covers a real fraction of the smaller segment;
    # a 1-bp nominal overlap is geometrically ADJACENT, not overlapping (v9.7.22-m).
    if overlap > 0 and ofrac >= OVERLAP_FRACTION_MIN:
        return "OVERLAPPING_REFERENCE_SEGMENTS", 0, overlap, ofrac, "COORDINATE", None
    gap = max(0, max(a0, b0) - min(a1, b1))
    if gap <= 5_000:
        return "ADJACENT_OR_NEARBY_REFERENCE_SEGMENTS", gap, overlap, ofrac, "COORDINATE", None
    if gap <= 25_000:
        return "DISTANT_ON_REFERENCE_CAUTION", gap, overlap, ofrac, "COORDINATE", None
    return "SAME_REFERENCE_ONLY_DISTANT", gap, overlap, ofrac, "COORDINATE", None


def _pair_score(a: ClusterBlastReference, b: ClusterBlastReference, bgc_a: BGCRecord, bgc_b: BGCRecord) -> tuple[int, str, dict[str, Any]]:
    # Component weights below are engineering judgment (not calibrated against a labelled
    # set). Precision is enforced by the downstream gate cascade (P-2, R1-R3, hub-degree),
    # not by these weights; the raw score is a coarse routing prior. HIGH>=14, MODERATE>=9.
    adjacency_class, gap, overlap, overlap_fraction, adjacency_basis, adjacency_span = _adjacency(a, b)
    score = 0
    if adjacency_class == "OVERLAPPING_REFERENCE_SEGMENTS":
        score += 8
    elif adjacency_class == "ADJACENT_OR_NEARBY_REFERENCE_SEGMENTS":
        score += 7
    elif adjacency_class == "DISTANT_ON_REFERENCE_CAUTION":
        score += 3
    elif adjacency_class == "SHARED_REFERENCE_NO_INTERVAL":
        score += 2
    else:
        score += 1

    if a.nprot and b.nprot:
        score += 3 if min(a.nprot, b.nprot) >= 3 else 1
        if (a.nprot + b.nprot) >= 12:
            score += 2
    if a.rank <= 5 and b.rank <= 5:
        score += 3
    elif a.rank <= 10 and b.rank <= 10:
        score += 1
    if a.mean_identity and b.mean_identity:
        avg_min = min(a.mean_identity, b.mean_identity)
        if avg_min >= 65:
            score += 2
        elif avg_min >= 50:
            score += 1
    else:
        avg_min = None

    type_overlap = sorted(_norm_type_tokens(a.reference_type) & _norm_type_tokens(b.reference_type))
    product_overlap = sorted(_product_tokens(bgc_a.products) & _product_tokens(bgc_b.products))
    if type_overlap:
        score += min(3, len(type_overlap))
    if product_overlap:
        score += min(2, len(product_overlap))
    if bgc_a.edge_status != "Interior" or bgc_b.edge_status != "Interior":
        score += 2
    if bgc_a.contig != bgc_b.contig:
        score += 1

    if score >= RGGMCI_SCORE_HIGH:
        conf = "HIGH_RG_GMCI_RESCUE"
    elif score >= RGGMCI_SCORE_MODERATE:
        conf = "MODERATE_RG_GMCI_CANDIDATE"
    else:
        conf = "LOW_SHARED_REFERENCE_SIGNAL"

    tiling = _subject_tiling(a.subjects, b.subjects)

    details = {
        "adjacency_class": adjacency_class,
        "adjacency_basis": adjacency_basis,
        "adjacency_gap": gap,
        "adjacency_span": adjacency_span,
        "gap_locus_suffix": gap,
        "overlap_locus_suffix": overlap,
        "overlap_fraction": overlap_fraction,
        "avg_min_identity": avg_min,
        "shared_reference_type_tokens": "; ".join(type_overlap),
        "shared_product_tokens": "; ".join(product_overlap),
        "subject_tiling_class": tiling["subject_tiling_class"],
        "n_shared_subjects": tiling["n_shared_subjects"],
        "n_a_only_subjects": tiling["n_a_only_subjects"],
        "n_b_only_subjects": tiling["n_b_only_subjects"],
        "a_only_subjects": tiling["a_only_subjects"],
        "b_only_subjects": tiling["b_only_subjects"],
        "shared_subjects": tiling["shared_subjects"],
    }
    return score, conf, details


def _geometry_gate(confidence: str, good_geometry_references: int, split_signature: bool) -> tuple[str, str]:
    """Geometry gates ACCEPTANCE, not just annotation (v9.7.22-m).

    good_geometry_references counts shared references whose segments OVERLAP or are ADJACENT on the
    reference. A pair can otherwise reach HIGH/MODERATE purely on convergence count (many shared
    references) with zero real geometry — the over-trust-on-convergence failure the v9.2 review named
    as the top structural hazard. Here a geometry-less pair is demoted: HIGH->MODERATE, and
    MODERATE->LOW unless a cross-scaffold physical-split signal independently supports it.
    """
    gate = "OK"
    if confidence == "HIGH_RG_GMCI_RESCUE" and good_geometry_references < 2:
        # v9.7.41: require >=2 good-geometry (OVERLAPPING/ADJACENT) references for HIGH. The 5-strain
        # PRIVATE calibration found 50-76% of HIGH pairs on every strain rode on a SINGLE good-geometry
        # reference -- the dominant over-promotion. A lone coarse-proxy adjacency is not enough for HIGH.
        confidence = "MODERATE_RG_GMCI_CANDIDATE"
        gate = "DEMOTED_HIGH_TO_MODERATE_fewer_than_2_good_geometry_references"
    if confidence == "MODERATE_RG_GMCI_CANDIDATE" and good_geometry_references == 0 and not split_signature:
        confidence = "LOW_SHARED_REFERENCE_SIGNAL"
        gate = "DEMOTED_TO_LOW_shared_reference_only_no_geometry_no_split_signal"
    return confidence, gate


def _same_contig(bgc_a: BGCRecord, bgc_b: BGCRecord) -> bool:
    return bool(bgc_a.contig) and contig_key(bgc_a.contig) == contig_key(bgc_b.contig)


def compute_rggmci(bgcs: list[BGCRecord], reference_map: dict[str, Any],
                   contigs: dict[str, str] | None = None,
                   glycosyltransferase_bgcs: set[str] | None = None,
                   contig_edge_bgcs: set[str] | None = None,
                   same_contig_pairs: bool = False) -> dict[str, Any]:
    # Only regions at a contig edge are paired; interior regions are whole on their contig (Alex, 2026-09-27).
    # Two regions on one contig are never a rescue: no assembly break separates them (Alex, 2026-09-27).
    edge_ids = _edge_status_fallback(bgcs) if contig_edge_bgcs is None else set(contig_edge_bgcs)
    all_records = [ClusterBlastReference(**r) for r in reference_map.get("reference_records", [])]
    # P-CBDB v9.7.100: subclusterblast finds sub-operon (cassette/sugar-operon) hits, NOT whole-cluster
    # hits, so it must NOT contribute to cluster-level rescue geometry. It is retained separately as a
    # functional marker (see subcluster_markers below) but excluded from the pairing pool.
    ref_records = [r for r in all_records if r.db_kind != "subclusterblast"]
    sub_records = [r for r in all_records if r.db_kind == "subclusterblast"]
    # v9.7.166 (#3): surface subclusterblast sub-operon hits per BGC as functional markers
    # (sugar-biosynthesis operons, PKS starters, NRPS monomers). Excluded from pairing geometry
    # above, but no longer dropped — attached to each BGC, top-5 by rank.
    _sub_by_bgc: dict[str, list] = {}
    for r in sub_records:
        _sub_by_bgc.setdefault(r.bgc_id, []).append(r)
    for _bid, _recs in _sub_by_bgc.items():
        _b = {b.bgc_id: b for b in bgcs}.get(_bid)
        if _b is not None and not _b.subcluster_hits:
            _b.subcluster_hits = [
                {"ref": rr.ref, "rank": rr.rank, "nprot": rr.nprot, "mean_identity": rr.mean_identity}
                for rr in sorted(_recs, key=lambda x: (x.rank if x.rank is not None else 999))[:5]
            ]
    by_bgc: dict[str, BGCRecord] = {b.bgc_id: b for b in bgcs}
    # Key the pool on (db_kind, ref): clusterblast and knownclusterblast hits to the same accession are
    # DIFFERENT evidence (cross-genome neighbor vs characterized MIBiG cluster) and are kept as separate
    # supporting references so a pair's evidence base records which DB types corroborate it.
    by_ref: dict[tuple[str, str], list[ClusterBlastReference]] = {}
    for rec in ref_records:
        by_ref.setdefault((rec.db_kind, rec.ref), []).append(rec)

    evidence_rows: list[dict[str, Any]] = []
    pair_acc: dict[tuple[str, str], dict[str, Any]] = {}

    # ── Phase 1: Evidence accumulation ──────────────────────────────────────────
    # For each shared reference, score every BGC pair that co-occurs on it.
    # Accumulate per-pair statistics across all shared references.
    for (db_kind, ref), recs in by_ref.items():
        # Skip references that only point to one query BGC.
        bgc_ids = sorted({r.bgc_id for r in recs if r.bgc_id in by_bgc and r.bgc_id in edge_ids})
        if len(bgc_ids) < 2:
            continue
        # Best record per BGC for this reference: rank first, then protein hits/score.
        best_by_bgc: dict[str, ClusterBlastReference] = {}
        for bgc_id in bgc_ids:
            candidates = [r for r in recs if r.bgc_id == bgc_id]
            best_by_bgc[bgc_id] = sorted(
                candidates,
                key=lambda r: (r.rank, -(r.nprot or 0), -(r.cumulative_score or 0)),
            )[0]
        # v9.7.409 (DEEP_AUDIT2_resource_dos #3): bound the within-strain all-vs-all BGC pairing.
        from .pair_scan_caps import cap_pair_scan_items
        for a_id, b_id in itertools.combinations(
                cap_pair_scan_items(sorted(best_by_bgc), label="rggmci"), 2):
            a = best_by_bgc[a_id]
            b = best_by_bgc[b_id]
            bgc_a, bgc_b = by_bgc[a_id], by_bgc[b_id]
            if not same_contig_pairs and _same_contig(bgc_a, bgc_b):
                continue
            score, confidence, details = _pair_score(a, b, bgc_a, bgc_b)
            row = {
                "pair": f"{a_id}+{b_id}",
                "bgc_a": a_id,
                "bgc_b": b_id,
                "ref": ref,
                "db_kind": db_kind,
                "source": a.source or b.source,
                "reference_type": a.reference_type or b.reference_type,
                "rank_a": a.rank,
                "rank_b": b.rank,
                "nprot_a": a.nprot,
                "nprot_b": b.nprot,
                "mean_identity_a": a.mean_identity,
                "mean_identity_b": b.mean_identity,
                "interval_a": f"{a.interval_start}-{a.interval_end}" if a.interval_start is not None and a.interval_end is not None else "",
                "interval_b": f"{b.interval_start}-{b.interval_end}" if b.interval_start is not None and b.interval_end is not None else "",
                **details,
            }
            if len(evidence_rows) < RGGMCI_MAX_EVIDENCE_ROWS:
                evidence_rows.append(row)

            key = tuple(sorted((a_id, b_id)))
            acc = pair_acc.setdefault(key, {
                "pair": f"{key[0]}+{key[1]}",
                "bgc_a": key[0],
                "bgc_b": key[1],
                "supporting_references": 0,
                "strong_supporting_references": 0,
                "complete_or_chromosome_references": 0,
                "good_geometry_references": 0,
                "mibig_good_geometry_references": 0,
                "scores": [],
                "sources": [],
                "reference_type_tokens": set(),
                "product_tokens": set(),
                "max_protein_sum": 0,
                "avg_min_identity_values": [],
                "complementary_disjoint_refs": 0,
                "overlapping_subject_refs": 0,
                "agg_a_only_subjects": set(),
                "agg_b_only_subjects": set(),
                "agg_shared_subjects": set(),
                "db_kinds": set(),
                "knownclusterblast_refs": 0,
                "clusterblast_refs": 0,
            })
            acc["supporting_references"] += 1
            acc["db_kinds"].add(db_kind)
            if db_kind == "knownclusterblast":
                acc["knownclusterblast_refs"] += 1
            elif db_kind == "clusterblast":
                acc["clusterblast_refs"] += 1
            if confidence == "HIGH_RG_GMCI_RESCUE":
                acc["strong_supporting_references"] += 1
            if "chromosome" in row["source"].lower() or row["ref"].startswith(("NZ_CP", "NC_", "CP")):
                acc["complete_or_chromosome_references"] += 1
            if row["adjacency_class"] in {"OVERLAPPING_REFERENCE_SEGMENTS", "ADJACENT_OR_NEARBY_REFERENCE_SEGMENTS"}:
                acc["good_geometry_references"] += 1
                # v9.7.42: of those, how many are characterized MIBiG clusters (BGC#######)? Cluster-level
                # co-location on a named MIBiG cluster is the corroboration the co-cluster exemption needs.
                if re.match(r"BGC\d{6,7}", row["ref"] or ""):
                    acc["mibig_good_geometry_references"] += 1
            acc["scores"].append(score)
            prot_sum = (a.nprot or 0) + (b.nprot or 0)
            acc["max_protein_sum"] = max(acc["max_protein_sum"], prot_sum)
            if details.get("avg_min_identity") is not None:
                acc["avg_min_identity_values"].append(details["avg_min_identity"])
            # Subject-gene tiling (v9.7.100). The evidence row's a_only/b_only are oriented to (a_id, b_id)
            # as iterated; the pair key is sorted, so re-orient onto (key[0], key[1]) before aggregating.
            tcls = details.get("subject_tiling_class")
            if tcls == "COMPLEMENTARY_DISJOINT":
                acc["complementary_disjoint_refs"] += 1
            elif tcls == "OVERLAPPING_SUBJECTS":
                acc["overlapping_subject_refs"] += 1
            if a_id == key[0]:
                acc["agg_a_only_subjects"].update(details.get("a_only_subjects", ()))
                acc["agg_b_only_subjects"].update(details.get("b_only_subjects", ()))
            else:
                acc["agg_a_only_subjects"].update(details.get("b_only_subjects", ()))
                acc["agg_b_only_subjects"].update(details.get("a_only_subjects", ()))
            acc["agg_shared_subjects"].update(details.get("shared_subjects", ()))
            if details.get("shared_reference_type_tokens"):
                acc["reference_type_tokens"].update(str(details["shared_reference_type_tokens"]).split("; "))
            if details.get("shared_product_tokens"):
                acc["product_tokens"].update(str(details["shared_product_tokens"]).split("; "))
            acc["sources"].append(
                f"{ref} ({row['adjacency_class']}, ranks {a.rank}/{b.rank}, proteins {(a.nprot or 0)}+{(b.nprot or 0)})"
            )

    # ── Phase 2: Pair scoring + initial confidence assignment ─────────────────
    # Aggregate per-pair evidence into a score, assign initial HIGH/MODERATE/LOW.
    ranked: list[dict[str, Any]] = []
    for key, acc in pair_acc.items():
        bgc_a = by_bgc[acc["bgc_a"]]
        bgc_b = by_bgc[acc["bgc_b"]]
        score = int(max(acc["scores"] or [0]) + min(6, acc["supporting_references"]) + min(4, acc["good_geometry_references"]) + min(3, acc["complete_or_chromosome_references"]))
        if score >= RGGMCI_SCORE_HIGH:
            confidence = "HIGH_RG_GMCI_RESCUE"
        elif score >= RGGMCI_SCORE_MODERATE:
            confidence = "MODERATE_RG_GMCI_CANDIDATE"
        else:
            confidence = "LOW_SHARED_REFERENCE_SIGNAL"
        split_signature = (bgc_a.edge_status == "Edge" and bgc_b.edge_status == "Edge"
                           and bool(bgc_a.contig) and bool(bgc_b.contig)
                           and contig_key(bgc_a.contig) != contig_key(bgc_b.contig))
        confidence, acceptance_gate = _geometry_gate(confidence, acc["good_geometry_references"], split_signature)
        # ── Phase 3: Sequential confidence gates (ORDER-DEPENDENT) ─────────────
        # Each gate demotes HIGH→MODERATE on a specific failure condition.
        # P-2 must run before co-cluster exemption (which un-demotes).
        # v9.7.41 P-2: product-class gate (demote-only) runs after geometry.
        conf_after_geometry = confidence
        confidence, product_gate = _product_class_gate(
            confidence, "; ".join(bgc_a.products), "; ".join(bgc_b.products))
        # v9.7.42: co-cluster carve-out — restore a complementary pair demoted only for incompatible
        # product classes when it is corroborated by >=2 ADJACENT/OVERLAPPING MIBiG cluster references.
        if confidence != conf_after_geometry and _cocluster_exempt(product_gate, acc["mibig_good_geometry_references"]):
            confidence = conf_after_geometry
            product_gate = "OK_cocluster_exempt_complementary_shared_mibig"
        if product_gate not in ("OK", "OK_shared_specific_product_class", "OK_compatible_hybrid",
                                "OK_cocluster_exempt_complementary_shared_mibig"):
            acceptance_gate = product_gate if acceptance_gate == "OK" else f"{acceptance_gate}+{product_gate}"
        elif acceptance_gate == "OK" and product_gate.startswith("OK_"):
            acceptance_gate = product_gate
        # R1 (v9.7.61): both-DROP-class pair suppression.
        confidence, r1_gate = _r1_drop_pair_gate(
            confidence, "; ".join(bgc_a.products), "; ".join(bgc_b.products))
        if r1_gate != "OK":
            acceptance_gate = r1_gate if acceptance_gate == "OK" else f"{acceptance_gate}+{r1_gate}"
        # R2 (v9.7.61): strong_supporting_references ≥ 1 co-requirement for HIGH.
        confidence, r2_gate = _r2_strong_ref_gate(confidence, acc["strong_supporting_references"])
        if r2_gate != "OK":
            acceptance_gate = r2_gate if acceptance_gate == "OK" else f"{acceptance_gate}+{r2_gate}"
        # R3 (v9.7.61): noise-class-only pair cap.
        confidence, r3_gate = _r3_noise_class_gate(
            confidence, "; ".join(bgc_a.products), "; ".join(bgc_b.products))
        if r3_gate != "OK":
            acceptance_gate = r3_gate if acceptance_gate == "OK" else f"{acceptance_gate}+{r3_gate}"
        avg_min = None
        if acc["avg_min_identity_values"]:
            avg_min = round(sum(acc["avg_min_identity_values"]) / len(acc["avg_min_identity_values"]), 1)
        # ── Phase 4: Subject-tiling + terminus-truncation rescue ──────────────
        # Pair-level subject-gene tiling verdict (v9.7.100, observational). Disjoint subject sets are the
        # DEFAULT for two unrelated BGCs, so a single disjoint reference is NOT evidence of a split — it must
        # be corroborated. COMPLEMENTARY_SPLIT requires (a) >= ST_MIN_DISJOINT_REFS shared references each
        # showing disjoint complementary tiling, (b) no shared reference showing the contigs co-hitting the
        # same subject genes, and (c) substantive tiling on both sides (>= ST_MIN_AGG_SUBJECTS_PER_SIDE
        # distinct subject genes aggregated per side). OVERLAPPING_PARALOG when shared subjects dominate
        # (same machinery, two loci). Otherwise INSUFFICIENT_SUBJECT_DATA. Observational only — no effect on
        # score/confidence in this cut.
        cd = acc["complementary_disjoint_refs"]
        ov = acc["overlapping_subject_refs"]
        n_a_side = len(acc["agg_a_only_subjects"])
        n_b_side = len(acc["agg_b_only_subjects"])
        substantive = (n_a_side >= ST_MIN_AGG_SUBJECTS_PER_SIDE and n_b_side >= ST_MIN_AGG_SUBJECTS_PER_SIDE)
        # Disjoint subjects are the DEFAULT for two unrelated BGCs, so the split verdict is only meaningful
        # where the geometry already supports linkage. Without >= ST_MIN_GEOMETRY_FOR_SPLIT good-geometry
        # references the pair is reported INSUFFICIENT regardless of disjointness — this is what stops the
        # verdict from labeling every low-signal shared-reference pair a "split". OVERLAPPING_PARALOG needs
        # no such gate: shared subject genes are positive evidence of paralogy on their own.
        geometry_ok = acc["good_geometry_references"] >= ST_MIN_GEOMETRY_FOR_SPLIT
        if ov >= 1 and cd == 0:
            subject_tiling_verdict = "OVERLAPPING_PARALOG"
        elif cd >= 1 and ov >= 1:
            subject_tiling_verdict = "MIXED_SUBJECT_SIGNAL"
        elif cd >= ST_MIN_DISJOINT_REFS and ov == 0 and substantive and geometry_ok:
            subject_tiling_verdict = "COMPLEMENTARY_SPLIT"
        else:
            subject_tiling_verdict = "INSUFFICIENT_SUBJECT_DATA"

        # P-TT v9.7.100: terminus-truncation rescue — physical, coordinate-first evidence that outranks
        # homology tiling. If either fragment is flush with its contig terminus, the two are on DIFFERENT
        # contigs, and they share a specific biosynthetic class, the cluster was sliced by the assembly
        # break. A short Full-contig partner (<= TT_SMALL_PARTNER_BP) that is entirely one BGC is the
        # classic broken-off arm. This OVERRIDES an OVERLAPPING_PARALOG verdict, because that verdict can be
        # an artifact of a gene class that legitimately occurs in multiple copies within ONE cluster
        # (e.g. bottromycin RRE/methyltransferase) appearing once on each severed fragment.
        tt_a = _terminus_truncation(bgc_a)
        tt_b = _terminus_truncation(bgc_b)
        tt_seq_a = _terminus_sequence_state(bgc_a, contigs)
        tt_seq_b = _terminus_sequence_state(bgc_b, contigs)
        low_complexity_junction = any(
            state in {"HOMOPOLYMER_TERMINUS", "LOW_COMPLEXITY_TERMINUS"}
            for state in (tt_seq_a, tt_seq_b)
        )
        diff_contig = bool(bgc_a.contig) and bool(bgc_b.contig) and contig_key(bgc_a.contig) != contig_key(bgc_b.contig)
        shared_class = _shared_class_tokens("; ".join(bgc_a.products), "; ".join(bgc_b.products))
        # Severed-arm geometry: a SHORT Full-contig fragment is a candidate broken-off arm. The rescue is
        # credible only when a terminus-truncated Edge region is paired with such a small complete contig
        # (the bottromycin NODE_69 case) OR the two share a SPECIFIC biosynthetic class (not a generic
        # umbrella). Either condition prevents the rule firing on every Edge region that happens to touch a
        # generic-token partner.
        a_small_full = bgc_a.edge_status == "Full-contig" and (bgc_a.contig_length or 1e9) <= TT_SMALL_PARTNER_BP
        b_small_full = bgc_b.edge_status == "Full-contig" and (bgc_b.contig_length or 1e9) <= TT_SMALL_PARTNER_BP
        # the truncated side must be the Edge fragment; the partner the small complete contig.
        severed_arm_geom = (
            (tt_a and b_small_full and bgc_a.edge_status == "Edge") or
            (tt_b and a_small_full and bgc_b.edge_status == "Edge")
        )
        terminus_truncation_rescue = bool(
            (tt_a or tt_b) and diff_contig and (bool(shared_class) or severed_arm_geom)
        )
        junction_evidence_state = "LONG_READ_ONLY" if low_complexity_junction else "SEQUENCE_COMPLEXITY_CLEAR"
        if terminus_truncation_rescue and low_complexity_junction:
            terminus_truncation_rescue = False
            override_note = "LONG_READ_ONLY_low_complexity_or_homopolymer_terminus"
        elif terminus_truncation_rescue:
            if subject_tiling_verdict == "OVERLAPPING_PARALOG":
                # Paralogy rests on a shared subject gene that can be multi-copy within one cluster; a
                # coordinate-confirmed truncation + severed-arm/specific-class geometry is stronger. Override.
                override_note = "OVERRODE_OVERLAPPING_PARALOG_terminus_truncation"
                subject_tiling_verdict = "TERMINUS_TRUNCATION_SPLIT"
            elif subject_tiling_verdict == "COMPLEMENTARY_SPLIT":
                override_note = "CORROBORATED_terminus_truncation"
            elif severed_arm_geom:
                # Only the strong severed-arm geometry (small complete partner contig) is allowed to PROMOTE
                # a weak/insufficient verdict to a split; a specific-class match alone annotates but does not
                # promote, to avoid the over-firing seen when generic geometry drove promotions.
                override_note = f"PROMOTED_FROM_{subject_tiling_verdict}_severed_arm"
                subject_tiling_verdict = "TERMINUS_TRUNCATION_SPLIT"
            else:
                override_note = "TERMINUS_FLAG_only_no_promotion"
        else:
            override_note = ""

        # ── Phase 4 gate (v9.7.352, AMBER_CORRECTNESS FIX 2): enforce the RG-GMCI two-proof rule ──
        # Governing principle (project memory `rggmci-two-proof`): a rescue needs homology geometry
        # (PROOF 1 — already gated by _geometry_gate/_product_class_gate/R1-R3 above) AND biosynthetic-
        # logic complementarity (PROOF 2 = subject_tiling_verdict). Until this cut the verdict was
        # OBSERVATIONAL and never gated confidence, so a pair the engine ITSELF labels OVERLAPPING_PARALOG
        # (same machinery on both contigs — affirmatively NOT a split, i.e. proof-2 REFUTED) still reached
        # HIGH and both BGCs collected the +8 triage rescue in scoring.py — the exact single-proof rescue
        # the two-proof rule exists to block. Because proof 2 is refuted (not merely absent), the pair is
        # demoted OUT of rescue eligibility entirely (below MODERATE) so it collects NO rescue bonus.
        # A terminus-truncation override has already re-labelled genuine severed arms to
        # TERMINUS_TRUNCATION_SPLIT above (coordinate-first physical evidence), so those are exempt.
        # COMPLEMENTARY_SPLIT / TERMINUS_TRUNCATION_SPLIT keep their confidence (proof 2 satisfied);
        # INSUFFICIENT_SUBJECT_DATA / MIXED_SUBJECT_SIGNAL are left untouched — absence or ambiguity of
        # subject data is not affirmative refutation, and demoting those would over-suppress the many
        # genuine data-sparse pairs the gate cascade already vets.
        if (subject_tiling_verdict == "OVERLAPPING_PARALOG"
                and confidence in ("HIGH_RG_GMCI_RESCUE", "MODERATE_RG_GMCI_CANDIDATE")):
            confidence = "LOW_SHARED_REFERENCE_SIGNAL"
            st_gate = "ST-PARALOG_no_complementarity_proof"
            acceptance_gate = st_gate if acceptance_gate == "OK" else f"{acceptance_gate}+{st_gate}"

        # Sugar-arm restore (see _sugar_arm_pair): only a pair the product-class gate alone demoted from HIGH, split
        # by its genes across contigs, tiling 2-RGGMCI_MAX_COCLUSTER_MIBIG_REFS characterized clusters together, with a
        # glycosyltransferase on the backbone side. R1-R3 still apply.
        if (confidence == "MODERATE_RG_GMCI_CANDIDATE" and conf_after_geometry == "HIGH_RG_GMCI_RESCUE"
                and product_gate == "DEMOTED_HIGH_TO_MODERATE_no_specific_product_class_after_exclusions"
                and diff_contig
                and subject_tiling_verdict in ("COMPLEMENTARY_SPLIT", "TERMINUS_TRUNCATION_SPLIT")
                and 2 <= acc["mibig_good_geometry_references"] <= RGGMCI_MAX_COCLUSTER_MIBIG_REFS
                and _sugar_arm_pair(bgc_a, bgc_b, glycosyltransferase_bgcs)):
            restored = "HIGH_RG_GMCI_RESCUE"
            pa, pb = "; ".join(bgc_a.products), "; ".join(bgc_b.products)
            restored, g1 = _r1_drop_pair_gate(restored, pa, pb)
            restored, g2 = _r2_strong_ref_gate(restored, acc["strong_supporting_references"])
            restored, g3 = _r3_noise_class_gate(restored, pa, pb)
            if restored == "HIGH_RG_GMCI_RESCUE":
                confidence = restored
                acceptance_gate = f"{acceptance_gate}+RESTORED_sugar_arm_glycosylated_partner"

        # P-CBDB v9.7.100: rescue evidence base — which DB types corroborate this pair. A rescue is most
        # credible when BOTH knownclusterblast (characterized-cluster identity) AND clusterblast (cross-genome
        # co-occurrence) support it; clusterblast-only still rescues (the AS-XXX/AS-XXX case where the BGC
        # has no characterized MIBiG hit but real genome neighbors); knownclusterblast-only is identity
        # without cross-genome corroboration. Observational; does not change score/confidence in this cut.
        dbk = acc["db_kinds"]
        if "knownclusterblast" in dbk and "clusterblast" in dbk:
            rescue_evidence_base = "BOTH_KCB_AND_CB"
        elif "knownclusterblast" in dbk:
            rescue_evidence_base = "KNOWNCLUSTERBLAST_ONLY"
        elif "clusterblast" in dbk:
            rescue_evidence_base = "CLUSTERBLAST_ONLY"
        else:
            rescue_evidence_base = "NONE"
        # split_signature + geometry acceptance gate computed above (v9.7.22-m).
        ranked.append({
            "pair": acc["pair"],
            "bgc_a": acc["bgc_a"],
            "contig_a": bgc_a.contig,
            "products_a": "; ".join(bgc_a.products),
            "edge_a": bgc_a.edge_status,
            "bgc_b": acc["bgc_b"],
            "contig_b": bgc_b.contig,
            "products_b": "; ".join(bgc_b.products),
            "edge_b": bgc_b.edge_status,
            "split_signature": split_signature,
            "acceptance_gate": acceptance_gate,
            "rggmci_score": score,
            "rggmci_confidence": confidence,
            "supporting_references": acc["supporting_references"],
            "strong_supporting_references": acc["strong_supporting_references"],
            "complete_or_chromosome_references": acc["complete_or_chromosome_references"],
            "good_geometry_references": acc["good_geometry_references"],
            "mibig_good_geometry_references": acc["mibig_good_geometry_references"],
            "avg_min_identity": avg_min,
            "rescue_evidence_base": rescue_evidence_base,
            "knownclusterblast_refs": acc["knownclusterblast_refs"],
            "clusterblast_refs": acc["clusterblast_refs"],
            "subject_tiling_verdict": subject_tiling_verdict,
            "terminus_truncation_rescue": terminus_truncation_rescue,
            "terminus_truncation_a": tt_a,
            "terminus_truncation_b": tt_b,
            "terminus_sequence_state_a": tt_seq_a,
            "terminus_sequence_state_b": tt_seq_b,
            "junction_evidence_state": junction_evidence_state,
            "terminus_override_note": override_note,
            "shared_class_tokens": "; ".join(sorted(shared_class)),
            "complementary_disjoint_refs": cd,
            "overlapping_subject_refs": ov,
            "n_a_only_subjects": len(acc["agg_a_only_subjects"]),
            "n_b_only_subjects": len(acc["agg_b_only_subjects"]),
            "n_shared_subjects": len(acc["agg_shared_subjects"]),
            "a_only_subjects": "; ".join(sorted(acc["agg_a_only_subjects"])[:25]),
            "b_only_subjects": "; ".join(sorted(acc["agg_b_only_subjects"])[:25]),
            "shared_subjects": "; ".join(sorted(acc["agg_shared_subjects"])[:25]),
            "max_protein_sum": acc["max_protein_sum"],
            "shared_reference_type_tokens": "; ".join(sorted(t for t in acc["reference_type_tokens"] if t)),
            "shared_product_tokens": "; ".join(sorted(t for t in acc["product_tokens"] if t)),
            "best_sources": " | ".join(acc["sources"][:5]),
            "interpretation_guard": "Homology-guided shared-reference linkage; not nucleotide-level joining.",
        })

    ranked.sort(key=lambda r: (r["rggmci_score"], r["supporting_references"], r["max_protein_sum"]), reverse=True)
    # ── Phase 5: Hub-degree guard + sort + split-candidate extraction ────────
    _apply_hub_degree_guard(ranked)
    high = [r for r in ranked if r["rggmci_confidence"] == "HIGH_RG_GMCI_RESCUE"]
    moderate = [r for r in ranked if r["rggmci_confidence"] == "MODERATE_RG_GMCI_CANDIDATE"]
    # split-signature candidates: cross-scaffold both-Edge pairs, ranked on their own so a genuine physical
    # split is not buried beneath higher-scoring same-contig / interior shared-reference pairs.
    split_candidates = sorted(
        [r for r in ranked if r["split_signature"]],
        key=lambda r: (r["rggmci_score"], r["supporting_references"], r["max_protein_sum"]),
        reverse=True,
    )
    basis_counts = {"COORDINATE": 0, "LOCUS_PROXY": 0, "NONE": 0}
    for r in evidence_rows:
        bk = r.get("adjacency_basis") or "NONE"
        basis_counts[bk] = basis_counts.get(bk, 0) + 1
    good_geom_pairs = sum(1 for r in ranked if r.get("good_geometry_references", 0) > 0)
    identity_rows = sum(1 for r in evidence_rows if r.get("avg_min_identity") is not None)
    complementary_split_pairs = sum(1 for r in ranked if r.get("subject_tiling_verdict") == "COMPLEMENTARY_SPLIT")
    overlapping_paralog_pairs = sum(1 for r in ranked if r.get("subject_tiling_verdict") == "OVERLAPPING_PARALOG")
    ceiling = "HIGH" if high else "MODERATE" if moderate else "LOW" if ranked else "NONE"
    summary_line = (
        f"RG-GMCI: {len(ranked)} pairs | good-geometry pairs {good_geom_pairs} | "
        f"complementary-split {complementary_split_pairs} / overlapping-paralog {overlapping_paralog_pairs} | "
        f"adjacency basis coord/proxy/none "
        f"{basis_counts['COORDINATE']}/{basis_counts['LOCUS_PROXY']}/{basis_counts['NONE']} | "
        f"identity populated {identity_rows}/{len(evidence_rows)} | promotion ceiling {ceiling}"
    )
    # Retain every HIGH-confidence pair even past the display cap, so a genuine split-pathway
    # rescue is never dropped by ranking truncation.
    _head = ranked[:RGGMCI_MAX_RANKED_PAIRS]
    _extra_high = [r for r in ranked[RGGMCI_MAX_RANKED_PAIRS:]
                   if r.get("rggmci_confidence") == "HIGH_RG_GMCI_RESCUE"]
    return {
        "status": "PASS" if ranked else "NULL_NO_RGGMCI_PAIRS",
        "algorithm": "RG-GMCI full pre-triage pass v1.9.5",
        "summary_line": summary_line,
        "rggmci_summary": {
            "pairs_total": len(ranked),
            "good_geometry_pairs": good_geom_pairs,
            "adjacency_basis_counts": basis_counts,
            "identity_rows_populated": identity_rows,
            "evidence_rows": len(evidence_rows),
            "complementary_split_pairs": complementary_split_pairs,
            "overlapping_paralog_pairs": overlapping_paralog_pairs,
            "promotion_ceiling": ceiling,
        },
        "pairing_scope": ("pairs_of_contig_edge_regions_on_different_contigs_sharing_reference_accession_from_clusterblast_txt"
                          if not same_contig_pairs else
                          "pairs_of_listed_regions_sharing_reference_accession_from_clusterblast_txt"),
        "contig_edge_regions": len(edge_ids & set(by_bgc)),
        "interior_regions_not_paired": len(set(by_bgc) - edge_ids),
        "reference_record_count": reference_map.get("reference_record_count", 0),
        "clusterblast_txt_file_count": len(reference_map.get("clusterblast_txt_files", [])),
        "pairs_total": len(ranked),
        "high_pairs": len(high),
        "moderate_pairs": len(moderate),
        "split_candidate_count": len(split_candidates),
        "split_candidates": split_candidates[:RGGMCI_MAX_RANKED_PAIRS],
        "ranked_pairs": _head + _extra_high,
        "evidence_rows": evidence_rows,
        "claim_safety": "RG-GMCI nominates split/linked BGC candidates from shared reference geometry. It must be completed before antibiotic/antifungal lead ranking, but it does not prove physical contig linkage or product identity.",
    }


RELATED_LOCUS_LABEL = "RELATED_LOCUS_NOT_A_RESCUE"


def related_locus_pairs(bgcs: list[BGCRecord], reference_map: dict[str, Any], contig_edge_bgcs: set[str],
                        contigs: dict[str, str] | None = None,
                        glycosyltransferase_bgcs: set[str] | None = None) -> list[dict[str, Any]]:
    """Pairs that involve an interior region or sit on one contig, kept apart from the rescue list (Alex, 2026-09-27).

    The same scoring runs over all regions; every HIGH or MODERATE pair with at least one region off a contig edge, or
    with both regions on the same contig, is returned with its grade moved to `shared_reference_grade` and
    `rggmci_confidence` set to RELATED_LOCUS_NOT_A_RESCUE. An interior region is whole on its contig, and two regions on
    one contig have no assembly break between them, so such a pair can point to related loci (a second copy, a separate
    locus of one pathway, a shared cassette) but never to a fragment the assembly separated. These pairs never feed
    the rescue confidence, candidate groups or triage. `related_because` says which rule placed the pair here.
    """
    by_id = {b.bgc_id: b for b in bgcs}
    everything = compute_rggmci(bgcs, reference_map, contigs=contigs, glycosyltransferase_bgcs=glycosyltransferase_bgcs,
                                contig_edge_bgcs=set(by_id), same_contig_pairs=True)
    out = []
    for p in everything.get("ranked_pairs", []):
        same = _same_contig(by_id[p["bgc_a"]], by_id[p["bgc_b"]])
        if p["bgc_a"] in contig_edge_bgcs and p["bgc_b"] in contig_edge_bgcs and not same:
            continue
        grade = p.get("rggmci_confidence", "")
        if not (grade.startswith("HIGH") or grade.startswith("MODERATE")):
            continue
        q = dict(p)
        q["shared_reference_grade"] = grade
        q["rggmci_confidence"] = RELATED_LOCUS_LABEL
        q["interior_side"] = "+".join(x for x, b in (("a", p["bgc_a"]), ("b", p["bgc_b"])) if b not in contig_edge_bgcs)
        q["related_because"] = "; ".join(r for r, hit in (("interior region", bool(q["interior_side"])),
                                                          ("same contig", same)) if hit)
        out.append(q)
    return out


def run_rggmci(zip_path: str | Path, bgcs: list[BGCRecord],
               contigs: dict[str, str] | None = None) -> dict[str, Any]:
    reference_map = parse_clusterblast_reference_map(zip_path, bgcs)
    gt = _glycosyltransferase_bgcs(zip_path, bgcs)
    edge = _contig_edge_bgcs(zip_path, bgcs)
    result = compute_rggmci(bgcs, reference_map, contigs=contigs, glycosyltransferase_bgcs=gt, contig_edge_bgcs=edge)
    result["related_locus_pairs"] = related_locus_pairs(bgcs, reference_map, edge, contigs=contigs,
                                                        glycosyltransferase_bgcs=gt)
    result["reference_map_status"] = reference_map.get("status")
    result["reference_parse_error_count"] = reference_map.get("parse_error_count", 0)
    return result
