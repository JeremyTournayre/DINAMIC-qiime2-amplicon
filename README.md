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