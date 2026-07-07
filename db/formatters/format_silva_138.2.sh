#!/usr/bin/env bash
# db/formatters/format_silva_138.2.sh
# Formats the SILVA 138.2 SSURef NR99 release into QIIME2 artifacts.
# Run this script once per QIIME2 environment before building a SILVA
# classifier under that environment (artifact formats can differ between
# QIIME2 releases, so ref-otus.qza/ref-taxonomy.qza must be rebuilt for
# each QIIME2 version you use).
#
# Author: Jeremy Tournayre
#
# Run from db/Data/ with the QIIME2 environment already active:
#   cd db/Data
#   source ../../config.sh && activate_qiime2
#   bash ../formatters/format_silva_138.2.sh
#
# Expected input files in SILVA_138.2_SSURef_NR99_tax_silva/:
#   tax_slv_ssu_138.2.txt              - taxonomy ranks
#   taxmap_slv_ssu_ref_nr_138.2.txt    - taxid map
#   tax_slv_ssu_138.2.tre              - taxonomy tree
#   SILVA_138.2_SSURef_NR99_tax_silva.fasta - reference sequences (RNA)
#
# Output: ref-otus.qza and ref-taxonomy.qza under a QIIME2-version-specific
# subdirectory, so re-running this formatter under a different QIIME2
# environment does not overwrite a previous version's artifacts:
#   SILVA_138.2_SSURef_NR99_tax_silva/<QIIME2_ENV>/ref-otus.qza
#   SILVA_138.2_SSURef_NR99_tax_silva/<QIIME2_ENV>/ref-taxonomy.qza
{

# CONDA_ENV_QIIME2 must already be set (config.sh sourced + activate_qiime2
# called by the person running this script — see usage note above).
if [[ -z "$CONDA_ENV_QIIME2" ]]; then
    echo "ERROR: CONDA_ENV_QIIME2 is not set. Run 'source ../../config.sh && activate_qiime2' first."
    exit 1
fi

export TMPDIR="tmp"
mkdir -p "$TMPDIR"

RAW="SILVA_138.2_SSURef_NR99_tax_silva"
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

echo "RUNNING: qiime import SILVATaxonomy"
qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path  $RAW/tax_slv_ssu_138.2.txt \
    --output-path $DB/taxranks_tax_slv_ssu_138.2.qza
echo "DONE: qiime import SILVATaxonomy"

echo "RUNNING: qiime import SILVATaxidMap"
qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path  $RAW/taxmap_slv_ssu_ref_nr_138.2.txt \
    --output-path $DB/taxmap_slv_ssu_ref_nr_138.2.qza
echo "DONE: qiime import SILVATaxidMap"

echo "RUNNING: qiime import Phylogeny"
qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path  $RAW/tax_slv_ssu_138.2.tre \
    --output-path $DB/taxtree_tax_slv_ssu_138.2.qza
echo "DONE: qiime import Phylogeny"

echo "RUNNING: qiime import RNASequence"
qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path  $RAW/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
    --output-path $DB/SILVA_138.2_SSURef_NR99_tax_silva.qza
echo "DONE: qiime import RNASequence"

echo "RUNNING: qiime rescript reverse-transcribe"
qiime rescript reverse-transcribe \
    --i-rna-sequences $DB/SILVA_138.2_SSURef_NR99_tax_silva.qza \
    --o-dna-sequences $DB/ref-otus.qza
echo "DONE: qiime rescript reverse-transcribe"

echo "RUNNING: qiime rescript parse-silva-taxonomy"
qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree   $DB/taxtree_tax_slv_ssu_138.2.qza \
    --i-taxonomy-map    $DB/taxmap_slv_ssu_ref_nr_138.2.qza \
    --i-taxonomy-ranks  $DB/taxranks_tax_slv_ssu_138.2.qza \
    --p-include-species-labels \
    --o-taxonomy        $DB/ref-taxonomy.qza
echo "DONE: qiime rescript parse-silva-taxonomy"

echo "Outputs: $DB/ref-otus.qza and $DB/ref-taxonomy.qza"
}