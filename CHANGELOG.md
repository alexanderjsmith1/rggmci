# Changelog

## 1.0.0rc2

- **Reference-guided completion.** For each region at a contig edge, the best KnownClusterBlast MIBiG reference
  (or, without one, a discovered one) is searched against every protein in the genome, not only those inside
  antiSMASH regions. It reports reference genes found elsewhere, genes split across two contig ends, and partner
  contigs with the tests each passed. It runs when a MIBiG protein database and DIAMOND or BLAST+ are found, and it
  never changes a score or a confidence.
- **`rggmci build-mibig-db`** builds that database from a folder of MIBiG GenBank files. It downloads nothing.
- New options: `--mibig-db`, `--reference-completion off`, `--sensitivity`, `--threads` and `--pfam`.
- Pairs gain `completion_tier`, `ref_completion_partner` and `split_gene_links`. With `--out-dir`, each genome also
  gets `_reference_completion.tsv`, `_split_genes.tsv` and `_partner_contigs.tsv`, and `SUMMARY.tsv` gains three
  completion columns.
- Without a database, or for a ZIP without its whole-genome GenBank file, the pairs are scored as before and the tier
  says why completion did not run.
- Pairs carry the read depth of both contigs, taken from SPAdes contig names, and `DEPTH_MISMATCH` when the lower is
  under 0.67 of the higher.
- **Residue tiling.** For pairs the paralog gate demoted, a DIAMOND search shows which stretch of each shared MIBiG
  protein each region covers. It runs with `--diamond-db`, or by default when the completion database has its DIAMOND
  index. It never changes a confidence.

## 1.0.0rc1

- First standalone release of RG-GMCI (see `src/rggmci/PROVENANCE.json` for the exact source version and
  the SHA-256 of every source file used).
- `rggmci` scores region pairs for one genome or a folder of genomes.
- `rggmci fasta` and `rggmci blastp-layer` add an optional BLASTp layer that you run yourself on NCBI.
- The result JSON records which antiSMASH ZIP and which build produced it.
