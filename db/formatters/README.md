# db/formatters/ — database-specific formatter scripts

Author: Jeremy Tournayre

Each script in this directory converts one reference database's native
format (SILVA, PR2, ...) into two QIIME2 artifacts:

```
db/Data/<Database>/<QIIME2_ENV>/ref-otus.qza
db/Data/<Database>/<QIIME2_ENV>/ref-taxonomy.qza
```

`<QIIME2_ENV>` is the QIIME2 conda environment name (e.g.
`qiime2-amplicon-2024.5`) that produced these artifacts. Segregating by
QIIME2 version is required: `.qza` files are QIIME2's own binary artifact
format, which can change between QIIME2 releases, so the same source
database formatted under two different QIIME2 versions produces two
genuinely different sets of files that must not overwrite one another.
`build_classifier.sh` expects this same layout and reads `ref-otus.qza` /
`ref-taxonomy.qza` from the subdirectory matching the currently active
`CONDA_ENV_QIIME2`.

This is different from analysis projects under `Data/` at the repository
root, which do not use a QIIME2-version subdirectory by default.

## How to run a formatter

All formatters use paths relative to `db/Data/`. Run them from that
directory, not from `db/` or the repository root. The QIIME2 environment
must be activated first: formatters read `$CONDA_ENV_QIIME2` to build their
output path, and they call `qiime tools import` directly — they do not
activate an environment themselves.

```bash
cd db/Data
source ../../config.sh
activate_qiime2
perl ../formatters/format_pr2_5.0.0.pl
bash ../formatters/format_pr2_5.0.0.sh
```

or, for a shell-only formatter:

```bash
cd db/Data
source ../../config.sh
activate_qiime2
bash ../formatters/format_silva_138.2.sh
```

This is a one-time manual step per database and per QIIME2 version you need
it under, since each database's raw download format differs and requires
source-specific parsing.

## Adding a new formatter

Copy the pattern of an existing formatter for a similar database:
- Keep all input paths relative to `db/Data/` (raw downloaded files stay at
  the top level; only the `qiime tools import` outputs go under
  `$CONDA_ENV_QIIME2/`).
- Read `$CONDA_ENV_QIIME2` from the environment rather than hardcoding a
  QIIME2 version or environment name.
- Do not hardcode absolute paths; set `TMPDIR` to a path under `db/Data/tmp`.
- Do not add a `conda activate` line inside the formatter itself.