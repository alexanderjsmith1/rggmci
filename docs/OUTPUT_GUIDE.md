# Reading rggmci output

## The input it needs

An antiSMASH result ZIP with its ClusterBlast output. RG-GMCI reads:

- the region GenBank files (`<contig>.region<NNN>.gbk`);
- the `knownclusterblast/` and `clusterblast/` text files. These are written when antiSMASH runs with
  KnownClusterBlast and ClusterBlast turned on (`--cb-knownclusters --cb-general` on the command line, or
  the matching options on the web server, where ClusterBlast is off by default).

Without the ClusterBlast files there is nothing to link, and every genome reports zero pairs. The
`reference_map_status` column says whether reference hits were found.

## Pairs (`<genome>_pairs.tsv`, or `--pairs`)

One row per scored pair of regions.

| column | meaning |
|---|---|
| `bgc_a`, `bgc_b` | region aliases, numbered in the order the regions were read |
| `contig_a`, `contig_b` | the full contig names antiSMASH used; use them to find the regions |
| `products_a`, `products_b` | antiSMASH's product annotations for each region |
| `edge_a`, `edge_b` | `Edge`: the region is within 5 kb of a contig end. `Full-contig`: it covers the contig. `Interior`: neither. `Unknown`: contig length not available |
| `rggmci_score` | the pair score before the acceptance gates |
| `rggmci_confidence` | `HIGH_RG_GMCI_RESCUE`, `MODERATE_RG_GMCI_CANDIDATE` or `LOW_SHARED_REFERENCE_SIGNAL`, after the gates |
| `acceptance_gate` | which gate kept or demoted the pair, and why |
| `supporting_references` | reference clusters hit by both regions |
| `good_geometry_references` | references where the two regions' hits overlap or sit side by side |
| `both_at_contig_ends` | whether both regions run to a contig end or cover their contig. A description of where the fragments sit, not a rating |
| `subject_tiling_verdict` | the gene-level check on shared references: `COMPLEMENTARY_SPLIT` (the regions hit different genes), `TERMINUS_TRUNCATION_SPLIT` (a fragment cut at a contig end), `OVERLAPPING_PARALOG` (they hit the same genes), `MIXED_SUBJECT_SIGNAL` (some of each) |
| `max_endpoint_hub_degree` | the larger number of HIGH partners either region had before the hub guard |
| `best_sources` | the five best supporting references, with where each region's hits sit on them and how many proteins each contributed |
| `depth_a`, `depth_b` | each contig's read depth, from the `_cov_` part of a SPAdes contig name; blank for other names |
| `depth_ratio`, `depth_flag` | the lower depth over the higher. `DEPTH_MISMATCH` below 0.67 (a plasmid, a repeat, or a second population), `DEPTH_CONSISTENT` otherwise, `DEPTH_UNAVAILABLE` without depths. A flag, never a demotion |
| `completion_tier` | whether reference-guided completion ran: `FULL`, `NO_ALIGNER`, `NO_MIBIG_PROTEINS`, `NO_WHOLE_GENOME_GENBANK` or `OFF` |
| `ref_completion_partner` | `yes` when one region's contig is an accepted partner for the other region's reference; `no`; `no_reference` when neither region had a reference. Blank unless `FULL` |
| `split_gene_links` | how many genes split across the two regions' contigs were called `CLEAR`; one broken gene counts once, however many references see it. Blank unless `FULL` |
| `residue_tiling` | for pairs the paralog gate demoted, which stretch of each shared MIBiG protein each region covers: `COMPLEMENTARY_RESIDUES`, `OVERLAPPING_RESIDUES`, `MIXED_RESIDUES`, `THIN_RESIDUES`, `NO_SHARED_MIBIG_REFERENCE` or `NOT_TESTED`. Present only when the DIAMOND search ran (see below) |

### `acceptance_gate` verdicts

A pair can carry several, joined by `+`.

| verdict | meaning |
|---|---|
| `OK_shared_specific_product_class` | the regions share a specific product class (compared by antiSMASH family) |
| `OK_compatible_hybrid` | the classes differ but form a known hybrid (e.g. NRPS with PKS) |
| `OK_cocluster_exempt_complementary_shared_mibig` | two different specific classes that would otherwise fail, kept because the regions tile 2 to 6 MIBiG clusters together |
| `RESTORED_sugar_arm_glycosylated_partner` | a sugar-only fragment returned to HIGH next to a backbone partner whose region names a glycosyltransferase |
| `DEMOTED_HIGH_TO_MODERATE_no_specific_product_class_after_exclusions` | no specific class left once saccharide, other, NAPAA, aminopolycarboxylic-acid and the RiPP umbrella word are set aside |
| `DEMOTED_HIGH_TO_MODERATE_incompatible_product_classes` | classes that do not combine into one pathway |
| `DEMOTED_HIGH_TO_MODERATE_distinct_ripp_subclasses` | two different RiPP subclasses |
| `DEMOTED_HIGH_TO_MODERATE_fewer_than_2_good_geometry_references` | fewer than two references where the hits overlap or sit side by side |
| `DEMOTED_HIGH_TO_MODERATE_R1_both_drop_class` | both regions carry only saccharide, NAPAA, hglE-KS, PKS-like or other |
| `DEMOTED_HIGH_TO_MODERATE_R2_no_strong_references` | no reference with strong support |
| `DEMOTED_HIGH_TO_MODERATE_R3_noise_class_both_sides` | neither region carries a specialist class (only fatty_acid, saccharide, NAPAA, hglE-KS or other) |
| `ST-PARALOG_no_complementarity_proof` | the gene-level check reads two copies, not two halves; the pair drops to LOW |
| `DEMOTED_HUB_PROMISCUITY_degree_N_gt_4` | one region has more than four HIGH partners; all its HIGH pairs step down |
| `DEMOTED_TO_LOW_shared_reference_only_no_geometry_no_split_signal` | the regions share references but nothing places them side by side |

**Only contig-edge regions are paired.** rggmci pairs a region only when antiSMASH flags it as touching a contig
edge (`/contig_edge="True"` on the region feature; the `Edge`/`Full-contig` status is used when a file has no flag).
A region antiSMASH places inside a contig is whole on that contig. No assembly break separated it, so rggmci cannot
override the assembly there. Two regions that both run to a contig end are fragments the assembly may have
separated; rggmci still does not claim a join.

## Related loci (`related_locus_pairs` in the JSON)

Pairs that involve an interior region, or two regions on one contig, are scored the same way and kept apart, with
`rggmci_confidence` set to `RELATED_LOCUS_NOT_A_RESCUE`, the original grade in `shared_reference_grade`, the interior
side in `interior_side`, and the reason in `related_because` (`interior region`, `same contig`, or both). They can
point to a second copy of a cluster, a separate locus of one pathway, or a shared cassette.
They are never contig rescues and never feed the rescue list, candidate groups or the hub count.

## Candidate groups (`CANDIDATE_GROUPS.tsv`)

HIGH pairs on different contigs that share a region are gathered into one group, so a pathway that may
be spread over several fragments reads as one candidate. `n_at_contig_ends` counts the group's regions that
run to a contig end. Groups are listed with the most HIGH pairs first.

MODERATE pairs never join groups. Chained together, they link a genome's similar clusters into one large
network that describes no single pathway. `possible_moderate_links` lists the cross-contig MODERATE pairs
that touch a group, as places to look next.

## Summary (`SUMMARY.tsv`)

One row per genome: regions, pairs, HIGH, MODERATE, how many HIGH pairs cross contigs and how many of those
have both regions at contig ends, how many groups, and the size of the largest group. Three columns describe
reference-guided completion: `completion_tier`, the number of accepted partner contigs, and the number of `CLEAR`
split genes. An `error` column records a ZIP that could not be read; the
batch carries on past it.

## Reference-guided completion

Written with `--out-dir` as three tables per genome, and under `reference_completion` in the JSON. Every row starts
with the tier, the aligner, the region (`<genome> / <contig> / regionNNN / <alias>`), its reference, how the reference
was chosen (`knownclusterblast`, or `discovered` with its protein count, ties and size check) and the reference's
compound names. When the tier is not `FULL`, each table holds one row that names the tier.

The reference is the region's best-ranked KnownClusterBlast MIBiG hit that the database holds. A region with no such
hit is given the MIBiG cluster with the most region proteins at 35% identity or more over half the protein or more:
at least two proteins, one of them with a biosynthetic role, in a cluster of 250 kb or less.

### Reference genes (`<genome>_reference_completion.tsv`)

One row per reference gene.

| column | meaning |
|---|---|
| `status` | `PRESENT_IN_CORE` (a match inside the region at 30% identity or more over half the reference protein, or 25% with e ≤ 1e-10), `MISSING_FOUND_CLEAR` (the best match elsewhere reaches 35% and beats the next candidate's score by 20% or more), `MISSING_FOUND_AMBIGUOUS`, `MISSING_NOT_FOUND` |
| `best_identity_pct`, `best_coverage_pct`, `best_locus`, `best_contig`, `best_region_identity` | the best match and where it sits; a match outside every region is labelled `(no antiSMASH region)` |
| `reciprocal_best` | whether this reference gene is that protein's own best match among the reference's genes |
| `partner_verdict` | for a find outside the region: `SUPPORTED` (another find of this reference sits within 10 kb), `SINGLE_GENE`, `HOUSEKEEPING_CONTEXT` (a lone find beside primary-metabolism genes; needs the Pfam scan) or `PARALOG_FAMILY` (the reciprocal search fails, or five or more other genome proteins match the same reference gene) |
| `reciprocal_best_mibig`, `reciprocal_best_identity`, `reciprocal_reference_identity`, `reciprocal_family_ratio` | the find's own best MIBiG cluster, and how its best score against the reference's compound family compares. Below 0.9 the find belongs to another family |
| `paralogs_in_genome`, `adjacent_finds`, `housekeeping_neighbours`, `neighbour_pfams` | the counts behind the verdict; `neighbour_pfams` lists only the housekeeping-marker domains found |

The reciprocal search uses the stretch of the genome protein that matched the reference and every MIBiG target, not
a top few. A conserved family can match many clusters almost equally.

### Split genes (`<genome>_split_genes.tsv`)

A reference gene in two pieces on different contigs. Each piece covers 40 residues or more of the gene and less than
80% of it. The pieces overlap by 20 residues or less and together cover half the gene or more. Each piece's open end
lies within 300 bp of a contig end.

| `split_call` | meaning |
|---|---|
| `CLEAR` | the pieces beat any whole-gene match elsewhere by 10 identity points or more |
| `WEAK` | they beat it by less |
| `RIVAL_STRONGER` | a whole gene elsewhere matches as well: likely paralog fragments |
| `MODULAR_UNRESOLVED` | the reference gene or a piece is an assembly-line protein (a KS or C domain, or two or more A domains), where module paralogy can fake two halves |
| `RECURRENT_COMMON_GENE` | the same two pieces also split against references of unrelated compound families, at least one belonging to a region that holds neither piece: a broken common gene, such as a regulator. `split_call_before_recurrence` keeps the call it replaced |

### Partner contigs (`<genome>_partner_contigs.tsv`)

One row per contig, other than the region's own, that carries a find or a split piece.

| column | meaning |
|---|---|
| `position_ok` | a find lies within 20 kb of a contig end, or the contig is under 40 kb |
| `depth_ratio`, `depth_ok` | the contig's read depth over the region's contig's; fails below 0.67, blank without depths (blank never fails) |
| `partner_verdicts` | the verdicts of its finds, and its split pieces |
| `support_ok` | two or more passing finds with one `SUPPORTED`, or a passing `CLEAR` split piece. A find passes when it is near a contig end and `SUPPORTED` or `SINGLE_GENE`; mobile-element genes (transposase, integrase, recombinase, insertion element) never count |
| `accepted`, `reason` | `accepted` needs support and a depth that does not fail. `reason` is `accepted`, `depth`, `support` or `no passing finds` |

### Residue tiling (pair fields beginning `residue_`)

For pairs the paralog gate demoted (`ST-PARALOG_no_complementarity_proof`), each region's proteins are aligned with
DIAMOND against the MIBiG clusters both regions share in KnownClusterBlast. Each side is placed on a cluster with
every residue used once. A cluster is complementary when each side places 300 residues or more and they share 10% or
less of the smaller side, and overlapping when they share more. It runs with `--diamond-db`, or by default when the
completion database has its DIAMOND index. With `--diamond-db`, `--residue-scope all` tests every pair. The detail fields describe the
cluster both sides match best. Against relatives near 50% identity, modular PKS often reads `MIXED_RESIDUES`: one
module cannot be told from another.

## BLASTp queries (`rggmci fasta`)

- `<genome>_blastp_NN.fasta`: one file per NCBI submission, at most 20 proteins and 85,000 residues.
  Header fields: genome, region alias, slot, role (`core`, `additional` or `edge`), gene, contig, region,
  coordinates, protein length and why it was chosen. Don't edit the headers; the layer reads them back.
- `blastp_manifest.json` and `blastp_manifest.tsv`: every query, its region and the pairs it serves, the
  selection settings, and the SHA-256 of the antiSMASH ZIP.
- `HOW_TO_BLASTP.md`: the steps on the NCBI site.

## BLASTp layer (`rggmci blastp-layer`)

- `blastp_proteins.tsv`, one row per query protein:
  - `blastp_status` is `HIT`, `NO_HITS` (XML2 lists the query with nothing), `NO_HITS_AFTER_EXCLUSION` or
    `NOT_IN_RESULTS` (not in any file given; could be not run or no hits);
  - the top hit's identity, query coverage, e-value, accession, title and organism;
  - `top_organisms` lists up to five, as named in the XML2.
- `blastp_pairs.tsv`, one row per pair:
  - the pair's `rggmci_confidence`, copied unchanged;
  - per side, proteins sent, with hits, with no hits, not in the results, and the median top identity;
  - `organisms_with_homologs_of_both` and `shared_organisms` need the XML2 and read "needs XML2" without it.
- `blastp_layer.json`: all of the above, plus query headers in the results that matched nothing in the
  manifest (renamed headers are reported, never guessed).

## What a candidate is, and what it is not

- RG-GMCI surfaces fragments by biosynthetic logic seen through homology. Two or more regions matched
  adjacent or overlapping parts of the same reference clusters, so they may belong to one pathway.
- That is one line of evidence. The second is the genes themselves: do the domains on the fragments
  complement each other into one coherent assembly line, with loading, extension, tailoring and release
  each accounted for once? Do they point to the same compound family?
- A candidate is not a joined contig. It does not show that the fragments are physically adjacent, and it
  does not name a compound.
- To check one, open each region's GenBank file and read the genes and domains. A long-read assembly
  settles physical adjacency.
