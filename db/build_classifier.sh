#!/usr/bin/env bash
# db/build_classifier.sh - QIIME2 naive-Bayes classifier builder.
#
# Author: Jeremy Tournayre
#
# Expects formatted reference sequences and taxonomy already in:
#   Data/<DATABASE>/<QIIME2_ENV>/ref-otus.qza
#   Data/<DATABASE>/<QIIME2_ENV>/ref-taxonomy.qza
# (produced by the database-specific formatter scripts in db/formatters/,
# run under the same QIIME2 environment).
#
# QIIME2_ENV is the active QIIME2 conda environment name (CONDA_ENV_QIIME2
# from config.sh). Classifiers are segregated by QIIME2 version because the
# same source database (e.g. SILVA) can be reformatted and reclassified
# under several QIIME2 releases, and .qza artifact formats can differ
# between releases.
#
# If primers are specified (PRIMER_F != "no"), extracts amplicon reads first.
# Output: Data/<DATABASE>/<QIIME2_ENV>/<DB_TARGET>/classifier.qza
#
# Usage (called by launch_db.pl):
#   bash build_classifier.sh <RUN_CONF>
{

# Load central configuration (CONDA_PROFILE, CONDA_ENV_QIIME2, etc.).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
activate_qiime2

source "$1"    # Load DB_TARGET, PRIMER_F, PRIMER_R, DATABASE, MIN_LENGTH, MAX_LENGTH.

export TMPDIR="tmp"

# All outputs for this database live under a QIIME2-version-specific
# subdirectory, since ref-otus.qza / ref-taxonomy.qza / classifier.qza are
# all produced by `qiime tools import` / `qiime feature-classifier`, and
# their exact binary format can change between QIIME2 releases.
DB_ROOT="Data/$DATABASE/$CONDA_ENV_QIIME2"
OUT_DIR="$DB_ROOT/$DB_TARGET"

echo "------------------------------------------"
echo "target       $DB_TARGET"
echo "primer_F     $PRIMER_F"
echo "primer_R     $PRIMER_R"
echo "database     $DATABASE"
echo "qiime2_env   $CONDA_ENV_QIIME2"
echo "min_length   $MIN_LENGTH"
echo "max_length   $MAX_LENGTH"

if [[ ! -f "$DB_ROOT/ref-otus.qza" ]] || [[ ! -f "$DB_ROOT/ref-taxonomy.qza" ]]; then
    echo "ERROR: $DB_ROOT/ref-otus.qza and/or ref-taxonomy.qza not found."
    echo "Run the appropriate formatter first (see db/formatters/README.md),"
    echo "with the SAME QIIME2 environment ($CONDA_ENV_QIIME2) activated."
    exit 1
fi

# Ensure the output directory exists even if this script is called directly
# (bypassing launch_db.pl, which normally creates it first for run.conf).
mkdir -p "$OUT_DIR"

# --- Record classifier provenance alongside this classifier's outputs --------
# Documents which QIIME2 environment, which pipeline commit, and which
# formatter run produced this classifier.
PROVENANCE_FILE="$OUT_DIR/PROVENANCE.txt"
{
    echo "Classifier build started: $(date -Iseconds)"
    echo "Database:      $DATABASE"
    echo "DB target:     $DB_TARGET"
    echo "QIIME2 env:    $CONDA_ENV_QIIME2"
    echo "Pipeline dir:  $PIPELINE_DIR"
    if git -C "$PIPELINE_DIR" rev-parse --is-inside-work-tree &>/dev/null; then
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

# If primers are provided, extract amplicon reads from the reference sequences.
ref_file="$DB_ROOT/ref-otus.qza"
if [[ $PRIMER_F != "no" ]]; then
    ref_file="$OUT_DIR/ref-seqs.qza"
    echo "RUNNING: feature-classifier extract-reads"
    qiime feature-classifier extract-reads \
        --i-sequences    "$DB_ROOT/ref-otus.qza" \
        --p-f-primer     $PRIMER_F \
        --p-r-primer     $PRIMER_R \
        --p-min-length   $MIN_LENGTH \
        --p-max-length   $MAX_LENGTH \
        --o-read-extraction-stats "$OUT_DIR/ref-seqs_extraction-stats.qza" \
        --o-reads        "$OUT_DIR/ref-seqs.qza"
    echo "DONE: feature-classifier extract-reads"
fi

echo "RUNNING: feature-classifier fit-classifier-naive-bayes"
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads     "$ref_file" \
    --i-reference-taxonomy  "$DB_ROOT/ref-taxonomy.qza" \
    --o-classifier          "$OUT_DIR/classifier.qza"
echo "DONE: feature-classifier fit-classifier-naive-bayes"

# Quick self-test: classify the training sequences.
echo "RUNNING: feature-classifier classify-sklearn (self-test)"
qiime feature-classifier classify-sklearn \
    --i-classifier      "$OUT_DIR/classifier.qza" \
    --i-reads           "$ref_file" \
    --o-classification  "$OUT_DIR/taxonomy.qza"
echo "DONE: feature-classifier classify-sklearn"

echo "RUNNING: metadata tabulate"
qiime metadata tabulate \
    --m-input-file      "$OUT_DIR/taxonomy.qza" \
    --o-visualization   "$OUT_DIR/taxonomy.qzv"
echo "DONE: metadata tabulate"

echo "------------DONE------------"
}