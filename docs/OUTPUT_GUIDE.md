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
have both regions at contig ends, how many groups, and the size of the largest group. An `error` column records a ZIP that could not be read; the
batch carries on past it.

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
