# rggmci: Reference-Guided Genome Mining Candidate Inference (RG-GMCI)

RG-GMCI surfaces antiSMASH regions in a fragmented genome assembly that may belong to one biosynthetic
pathway: fragments that the assembly left on separate contigs. It reads biosynthetic logic through the
ClusterBlast and KnownClusterBlast results that antiSMASH already wrote.

It is a homology-guided candidate method, not a contig joiner. If the reads had supported a join, the
assembler would have made it. A candidate is a hypothesis to check at gene level.

It has no required dependencies. Biopython is used when installed; otherwise a built-in
GenBank reader is used.

## Run it

```bash
pip install rggmci-*.whl
rggmci my_genome.antismash.zip --out result.json --pairs pairs.tsv     # one genome
rggmci antismash_zips/ --out-dir rggmci_results/                        # many genomes
```

A batch writes one JSON and one pairs table per genome, plus `SUMMARY.tsv` and `CANDIDATE_GROUPS.tsv`.
[docs/OUTPUT_GUIDE.md](docs/OUTPUT_GUIDE.md) explains every column.

The input is an antiSMASH result ZIP that includes its `knownclusterblast/` and `clusterblast/` output.
A ZIP with only region GenBank files still runs, but no pairs can be scored without the reference hits.

### Starting from an assembly

rggmci does not call genes or detect clusters. antiSMASH does both, so run it on the assembly first:
- **Web server** ([antismash.secondarymetabolites.org](https://antismash.secondarymetabolites.org)): upload the
  FASTA and turn on **ClusterBlast**, which is off by default (or click "All on"). KnownClusterBlast is on by
  default. Download the result ZIP.
- **Command line:** `antismash --genefinding-tool prodigal --cb-general --cb-knownclusters --output-dir out assembly.fasta`,
  then give rggmci the result ZIP antiSMASH writes in `out/`.

The proteins `rggmci fasta` writes are antiSMASH's own gene calls, so no separate gene calling or Pfam scan
is needed.

## A second layer from your own BLASTp

`rggmci fasta` writes the proteins of the regions in HIGH cross-contig pairs as FASTA files sized for NCBI
web BLASTp, plus a manifest and a short how-to. You run the searches yourself. `rggmci blastp-layer` reads the
results back.

```bash
rggmci fasta my_genome.antismash.zip --out-dir queries/          # then BLASTp each queries/*.fasta on NCBI
rggmci blastp-layer --manifest queries/blastp_manifest.json \
       --hits Hit_Table.csv --xml2 results.xml --out-dir layer/
```

- Per region it takes up to 6 biosynthetic genes from the whole region, core genes first and nearest the
  contig end first. `--edge-genes N` adds the N genes nearest the break, of any kind. `--include-moderate` adds
  MODERATE pairs.
- Download both the Hit Table (CSV) and the Single-file XML2. The XML2 carries organism names and lists the
  queries that found nothing. A Hit Table alone cannot tell "no hits" from "not run".
- The layer reports, per pair, which source organisms carry homologs of both fragments' proteins. If the
  genome is already in NCBI, pass `--exclude-organism "<its species>"`, or its own proteins come back as hits.
- The layer sits beside the RG-GMCI confidence and never changes it. Similarity is not identity: shared
  homologs in one organism are a reason to look at that organism's genome, not proof of one pathway.

## What it reports

- One row per scored region pair, with the two regions' contigs, products and edge status.
- `rggmci_confidence`: `HIGH_RG_GMCI_RESCUE`, `MODERATE_RG_GMCI_CANDIDATE` or
  `LOW_SHARED_REFERENCE_SIGNAL`.
  - HIGH needs a pair score over the threshold. The score counts shared references, references where
    the two regions' hits overlap or sit side by side, and whether the regions touch a contig end or
    lie on different contigs.
  - HIGH also needs at least two such overlapping or side-by-side references, and a specific shared
    product class or a known compatible hybrid.
  - Only regions antiSMASH flags as touching a contig edge are paired, and never two regions on one contig. A region
    antiSMASH places inside a contig is whole there, and two regions on one contig have no assembly break between
    them; such pairs are listed apart as related loci (`related_locus_pairs` in the JSON), never as rescues.
  - Product classes are compared by antiSMASH family, so an arylpolyene or PKS-like fragment matches a type I
    PKS as PKS, and a lanthipeptide matches as RiPP.
  - Further gates demote a pair that fails them: too few strong references, noise-only classes, a
    region linked to more than four HIGH partners, or paralog-like overlap. `acceptance_gate` says which gate
    acted; [docs/OUTPUT_GUIDE.md](docs/OUTPUT_GUIDE.md) lists every verdict.
  - Paralog-like overlap means the two regions hit the same genes of a reference. Up to 2 shared genes are
    tolerated when they are at most 15% of all genes hit and each region hits 2 genes of its own.
  - A pair held back only by the product-class gate, where one fragment carries only sugar genes, can return to
    HIGH when its partner is a backbone class (PKS,
    NRPS, RiPP, terpene) whose own region file names a glycosyltransferase. The two must be on different
    contigs, pass the gene-level split check, and tile 2 to 6 MIBiG clusters together.
- The references that support each pair.
- `both_at_contig_ends`: whether both regions of a pair run to a contig end. It describes where the
  fragments sit and hides nothing.
- Candidate groups: HIGH pairs on different contigs that share regions, gathered into one group. That is
  how a pathway spread over three or more fragments shows up.

Regions are named by contig and region number, exactly as antiSMASH named them, so every pair can be
found in the antiSMASH output. Cite a region by contig and region, never by its `BGC` alias alone: the alias is
rggmci's reading order and changes between antiSMASH runs.

[docs/WORKED_EXAMPLES.md](docs/WORKED_EXAMPLES.md) walks through four pairs from public genomes: a complementary
split, two copies of one cluster, a hub, and a pair with too little support.

## What it does not do

- It does not join contigs, fill gaps, or show that two fragments are physically adjacent.
- It does not name a compound. ClusterBlast similarity is similarity, not identity.
- Region products come from antiSMASH's own annotations. A shared product label is a reason to look,
  not proof of one pathway.

## Install

```bash
pip install "git+https://github.com/alexanderjsmith1/rggmci"          # or: pip install rggmci-*.whl
pip install "rggmci[bio] @ git+https://github.com/alexanderjsmith1/rggmci"  # with Biopython
```

Python 3.12 or newer. No other dependency is required.

## Cite

If you use rggmci, cite it (see `CITATION.cff`) and cite antiSMASH, whose output it reads:
antiSMASH 8.0, doi:10.1093/nar/gkaf334. Reference clusters come from MIBiG 4.0, doi:10.1093/nar/gkae1115.

## Where it comes from

`src/rggmci/PROVENANCE.json` records the source version this build came from, and the SHA-256 of every
source file used.

## Licence

MIT. See `LICENSE`.
