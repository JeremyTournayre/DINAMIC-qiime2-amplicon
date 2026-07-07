#!/usr/bin/env bash
# db/formatters/format_pr2_5.0.0.sh
# Imports the reformatted PR2 v5.0.0 files into QIIME2 artifacts.
# Run format_pr2_5.0.0.pl first to generate the input files.
#
# Author: Jeremy Tournayre
#
# Run from db/Data/ with the QIIME2 environment already active:
#   cd db/Data
#   source ../../config.sh && activate_qiime2
#   perl ../formatters/format_pr2_5.0.0.pl
#   bash ../formatters/format_pr2_5.0.0.sh
#
# Output: ref-otus.qza and ref-taxonomy.qza under a QIIME2-version-specific
# subdirectory, so re-running this formatter under a different QIIME2
# environment does not overwrite a previous version's artifacts:
#   pr2_version_5.0.0/<QIIME2_ENV>/ref-otus.qza
#   pr2_version_5.0.0/<QIIME2_ENV>/ref-taxonomy.qza
{

if [[ -z "$CONDA_ENV_QIIME2" ]]; then
    echo "ERROR: CONDA_ENV_QIIME2 is not set. Run 'source ../../config.sh && activate_qiime2' first."
    exit 1
fi

export TMPDIR="tmp"
mkdir -p "$TMPDIR"

RAW="pr2_version_5.0.0"
DB="$RAW/$CONDA_ENV_QIIME2"
mkdir -p "$DB"

# --- Record provenance for this formatted database --------------------------
# Documents which QIIME2 environment produced ref-otus.qza / ref-taxonomy.qza
# under this directory — the directory name alone (the QIIME2 env) is fragile
# if the folder is ever moved, renamed, or copied elsewhere; this file makes
# the origin explicit and durable.
PROVENANCE_FILE="$DB/PROVENANCE.txt"
{
    echo "Database formatting started: $(date -Iseconds)"
    echo "Database:     $RAW"
    echo "QIIME2 env:   $CONDA_ENV_QIIME2"
    echo "Formatter:    $(basename "${BASH_SOURCE[0]}")"
    if [[ -n "$PIPELINE_DIR" ]] && git -C "$PIPELINE_DIR" rev-parse --is-inside-work-tree &>/dev/null; then
        echo "Git commit:  $(git -C "$PIPELINE_DIR" rev-parse HEAD)"
        echo "Git tag:     $(git -C "$PIPELINE_DIR" describe --tags --always)"
        if [[ -n "$(git -C "$PIPELINE_DIR" status --porcelain)" ]]; then
            echo "WARNING: working directory has uncommitted changes at run time!"
            git -C "$PIPELINE_DIR" status --porcelain
        fi
    else
        echo "WARNING: pipeline directory is not a git repository (no provenance available)."
    fi
    echo "---"
} >> "$PROVENANCE_FILE"

echo "RUNNING: qiime import Sequence (PR2)"
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path  $RAW/formatted_pr2_version_5.0.0_SSU_UTAX.fasta \
    --output-path $DB/ref-otus.qza
echo "DONE: qiime import Sequence"

echo "RUNNING: qiime import Taxonomy (PR2)"
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path   $RAW/taxonomy_pr2_version_5.0.0_SSU_UTAX.fasta \
    --output-path  $DB/ref-taxonomy.qza
echo "DONE: qiime import Taxonomy"

echo "Outputs: $DB/ref-otus.qza and $DB/ref-taxonomy.qza"
}