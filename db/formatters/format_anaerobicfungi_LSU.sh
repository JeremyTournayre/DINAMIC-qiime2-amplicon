#!/usr/bin/env bash
# db/formatters/format_anaerobicfungi_LSU.sh
# Formats the AGF LSU (anaerobic gut fungi, LSU) database v2 into QIIME2 artifacts.
# Run this script once per QIIME2 environment before building an AGF LSU
# classifier under that environment (artifact formats can differ between
# QIIME2 releases, so ref-otus.qza/ref-taxonomy.qza must be rebuilt for
# each QIIME2 version you use).
#
# Author: Jeremy Tournayre
#
# Run from db/Data/ with the QIIME2 environment already active:
#   cd db/Data
#   source ../../config.sh && activate_qiime2
#   bash ../formatters/format_anaerobicfungi_LSU.sh
#
# Expected input files in anaerobicfungi_LSU/:
#   AGF_LSU_database_v2.fasta - reference sequences (DNA)
#   AGF_LSU_tax_v2.txt        - taxonomy (headerless TSV)
#
# Output: ref-otus.qza and ref-taxonomy.qza under a QIIME2-version-specific
# subdirectory, so re-running this formatter under a different QIIME2
# environment does not overwrite a previous version's artifacts:
#   anaerobicfungi_LSU/<QIIME2_ENV>/ref-otus.qza
#   anaerobicfungi_LSU/<QIIME2_ENV>/ref-taxonomy.qza
{

# CONDA_ENV_QIIME2 must already be set (config.sh sourced + activate_qiime2
# called by the person running this script — see usage note above).
if [[ -z "$CONDA_ENV_QIIME2" ]]; then
    echo "ERROR: CONDA_ENV_QIIME2 is not set. Run 'source ../../config.sh && activate_qiime2' first."
    exit 1
fi

export TMPDIR="tmp"
mkdir -p "$TMPDIR"

RAW="anaerobicfungi_LSU_v2"
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

# --- Clean the raw download before import ------------------------------------
# The v2 files have Windows line endings, and 4 fasta IDs differ from the
# taxonomy file (WI3B vs WI3_B, ITS158SITS2 vs ITS1_58S_ITS2, UH31 vs UH3_1).
# fit-classifier-naive-bayes silently drops reads that have no taxonomy, and
# classify-consensus-vsearch refuses to run, so fix the IDs and check them.
# The downloaded files in $RAW are left untouched.
CLEAN="$DB/clean"
mkdir -p "$CLEAN"
FASTA="$CLEAN/AGF_LSU_database_v2.fasta"
TAX="$CLEAN/AGF_LSU_tax_v2.txt"
tr -d '\r' < "$RAW/AGF_LSU_database_v2.fasta" | sed \
    -e 's/^>Neocallimastix_sp_WI3B_ITS158SITS2_/>Neocallimastix_sp_WI3_B_ITS1_58S_ITS2_/' \
    -e 's/^>Neocallimastix_sp_WI3B_LSU_/>Neocallimastix_sp_WI3_B_LSU_/' \
    -e 's/^>Piromyces_spUH31_LSU_/>Piromyces_spUH3_1_LSU_/' > "$FASTA"
tr -d '\r' < "$RAW/AGF_LSU_tax_v2.txt" | sed -e '$a\' > "$TAX"
ID_DIFF=$(comm -3 <(grep '>' "$FASTA" | sed 's/^>//; s/[[:space:]].*//' | sort) \
                  <(cut -f1 "$TAX" | sort))
if [[ -n "$ID_DIFF" ]]; then
    echo "ERROR: sequence IDs differ between fasta and taxonomy:"
    echo "$ID_DIFF"
    exit 1
fi
echo "IDs checked: $(grep -c '>' "$FASTA") sequences, all with a taxonomy" | tee -a "$PROVENANCE_FILE"

echo "RUNNING: qiime import Sequence"
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path  "$FASTA" \
    --output-path $DB/ref-otus.qza
echo "DONE: qiime import Sequence"

echo "RUNNING: qiime import Taxonomy"
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path  "$TAX" \
    --output-path $DB/ref-taxonomy.qza
echo "DONE: qiime import Taxonomy"

echo "Outputs: $DB/ref-otus.qza and $DB/ref-taxonomy.qza"
}
