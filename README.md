# Holoruminant QIIME2 amplicon pipeline — supplementary scripts

Author: Jeremy Tournayre

This repository contains the scripts used to process 16S/18S amplicon
sequencing data (QIIME2, DADA2, ANCOM-BC2, DESeq2, LEfSe, IndicSpecies).

## Requirements

- QIIME2 (amplicon distribution), installed via conda
- R with the following packages: `qiime2R`, `phyloseq`, `vegan`,
  `pairwiseAdonis`, `ape`, `ANCOMBC`, `DESeq2`, `indicspecies`, `lefser`,
  `cowplot`, `Cairo`, `data.table`, `doSNOW`, `doRNG`, `microbiome`,
  `openxlsx`, `qiime2R`
- Perl with `Text::CSV`
- `dos2unix`, `unzip`

All machine-specific paths (conda environment name, R library path,
database location) are set once in `config.sh`.

## Repository structure

```
config.sh              Central configuration (paths, QIIME2 environment)
launch.pl               Orchestrator: reads input.txt, runs qiime2.sh per target
launch_multipleqsub.pl  Submits one launch.pl job per target on an SGE cluster
qiime2.sh               Main per-target pipeline (import through diversity/DA)
bin/R/                  R helper scripts (rarefaction depth, diversity, plots)
bin/perl/               Perl helper scripts (sample filtering, result merging)
db/                     Reference database formatting and classifier building
db/formatters/          One script per reference database (SILVA, PR2)
```

## Usage

### Running the main pipeline

1. Edit `config.sh` to match your environment (QIIME2 conda environment
   name, R library path).
2. Prepare `input.txt` (tab-separated, one row per analysis) and, for each
   target, a `Data/<target>/` folder containing `input/` (raw fastq.gz) and
   `mapping.tsv` (sample metadata, QIIME2-compatible types).
3. Run:
   ```bash
   perl launch.pl input.txt
   ```
   or, to submit one job per row on a cluster:
   ```bash
   source config.sh && activate_qiime2
   perl launch_multipleqsub.pl input.txt
   ```

Targets are created under `Data/<qiime2_conda_env>/<target>/`. Each run
writes a `PROVENANCE.txt` file recording the exact pipeline git commit used,
for traceability.

### vsearch chimera thresholds (optional columns)

When `vsearch` is set to `vsearch`, `qiime vsearch uchime-denovo` runs after
DADA2. Its thresholds are read from the optional `input.tsv` columns
`vsearch_minh` and `vsearch_mindiv`; when absent or empty, the QIIME2 defaults
are used (0.28 / 0.8). Runs before 2026-11-13 used 0.05 / 0.2 hard-coded
(values from a fungal LSU study, too strict for bacterial 16S: 20-45 % of
DADA2 reads removed); set these values to reproduce older results.

### Blocking factor and statistics (optional `block` column)

The optional `input.tsv` column `block` names a `mapping.tsv` column (e.g.
`Animal` for the donor animal) used as a covariate in the statistics. Empty
or absent: no block. It is ignored for an analysis whose comparison column
is the block itself, or when block and group are confounded.

With a block:

- QIIME2: `qiime diversity adonis --p-formula "<block>+<comparison>"`
  (free permutations), `adonis_<block>_<metric>.qzv`.
- vegan (`diversity_vegan_TSS/`): `permanova_block_*.tsv`
  (`adonis2(dist ~ block + group, by = "margin")`, free and within-block
  permutations), `pairwise_block_*.tsv` (> 2 groups), `alpha_lm_*.tsv`
  (`alpha ~ block + group`, marginal F tests).
- ANCOM-BC2 `fix_formula = "<block> + <group>"`, DESeq2
  `~ <block> + <group>`, IndicSpecies permutations within blocks. LEfSe has
  no block.

Always (since 2026-11):

- Aitchison distance (CLR, pseudocount 1, non-rarefied table) in QIIME2
  (`aitchison_distance_matrix.qza/.qzv`, PCoA) and vegan; dispersion tests
  (`permdisp_<metric>.qzv`, `betadisper_*.tsv`).
- ANCOM-BC2 `struc_zero = TRUE` (no manual pseudocount); `neg_lb` is FALSE
  unless every group has > 30 samples. Structural zeros are listed in
  `ANCOM-BC_structural_zeros.tsv`.
- IndicSpecies and LEfSe on relative abundances; IndicSpecies `p.value.BH`.
- DESeq2 pairwise log fold changes are log2(group2 / group1).

### Additional numeric metadata (mapping_extra.tsv)

Some numeric metadata columns (e.g. physico-chemical measurements) are
sparse enough (many missing values) to make QIIME2's own metadata-driven
commands fail. These columns are kept in a separate file,
`Data/<target>/mapping_extra.tsv` (same format as `mapping.tsv`, indexed by
the same sample IDs), which is only read by `add_metadata_on_diversity_plot.R`
to project numeric variables onto ordination plots. `mapping.tsv` itself
stays limited to the columns QIIME2 needs.

### Building a reference database classifier

1. Download and format the raw database with the matching script in
   `db/formatters/` (see `db/formatters/README.md`).
2. Prepare `db/input.tsv` (one row per primer set) and run:
   ```bash
   cd db
   source config.sh && activate_qiime2
   perl launch_db.pl input.tsv
   ```

Classifiers are stored under
`db/Data/<Database>/<qiime2_conda_env>/<Primer>/classifier.qza`.

## Notes

- The pipeline is resumable: each step is skipped if its output already
  exists. If you change an input parameter or file (e.g. a filter
  definition, `mapping.tsv`), delete the corresponding output file(s)
  before rerunning so they are regenerated.
- Random seeds are fixed (`set.seed(42)`) in R scripts using permutation-
  based tests (PERMANOVA, IndicSpecies) for reproducibility.