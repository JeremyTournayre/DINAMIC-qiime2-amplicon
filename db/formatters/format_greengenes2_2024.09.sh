#!/usr/bin/env bash
# db/formatters/format_greengenes2_2024.09.sh
# Formats the Greengenes2 2024.09 backbone (full-length 16S) into QIIME2
# artifacts. Run this script once per QIIME2 environment before building a
# Greengenes2 classifier under that environment.
#
# Author: Jeremy Tournayre
#
# Run from db/Data/ with the QIIME2 environment already active:
#   cd db/Data
#   source ../../config.sh && activate_qiime2
#   bash ../formatters/format_greengenes2_2024.09.sh
#
# The backbone is used rather than the full ASV set (2024.09.seqs.fna.qza):
# training on the full set takes more than a week for near-identical results
# (backbone: ~1 day).
#
# Expected input files in GreenGenes2_2024_09_backbone/ (downloaded if absent):
#   2024.09.backbone.full-length.fna.qza - backbone sequences
#   2024.09.backbone.tax.qza             - backbone taxonomy
# Source: https://ftp.microbio.me/greengenes_release/2024.09/
#
# These .qza are distributed by Greengenes2 itself (built with an older QIIME2
# release), so they are exported and re-imported here under the active
# QIIME2 environment, like the other formatters.
#
# Output: ref-otus.qza and ref-taxonomy.qza under a QIIME2-version-specific
# subdirectory:
#   GreenGenes2_2024_09_backbone/<QIIME2_ENV>/ref-otus.qza
#   GreenGenes2_2024_09_backbone/<QIIME2_ENV>/ref-taxonomy.qza
{

# CONDA_ENV_QIIME2 must already be set (config.sh sourced + activate_qiime2
# called by the person running this script — see usage note above).
if [[ -z "$CONDA_ENV_QIIME2" ]]; then
    echo "ERROR: CONDA_ENV_QIIME2 is not set. Run 'source ../../config.sh && activate_qiime2' first."
    exit 1
fi

set -e

export TMPDIR="tmp"
mkdir -p "$TMPDIR"

RAW="GreenGenes2_2024_09_backbone"
DB="$RAW/$CONDA_ENV_QIIME2"
URL="https://ftp.microbio.me/greengenes_release/2024.09"
SEQS="2024.09.backbone.full-length.fna.qza"
TAX="2024.09.backbone.tax.qza"
mkdir -p "$DB"

# --- Record provenance for this formatted database --------------------------
PROVENANCE_FILE="$DB/PROVENANCE.txt"
{
    echo "Database formatting started: $(date -Iseconds)"
    echo "Database:     $RAW"
    echo "Source:       $URL/$SEQS , $URL/$TAX"
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

echo "RUNNING: download Greengenes2 2024.09 backbone"
for f in "$SEQS" "$TAX"; do
    [[ -f "$RAW/$f" ]] || wget -O "$RAW/$f" "$URL/$f"
done
md5sum "$RAW/$SEQS" "$RAW/$TAX" >> "$PROVENANCE_FILE"
echo "DONE: download Greengenes2 2024.09 backbone"

echo "RUNNING: qiime tools export"
qiime tools export --input-path $RAW/$SEQS --output-path $RAW/export_seqs
qiime tools export --input-path $RAW/$TAX  --output-path $RAW/export_tax
echo "DONE: qiime tools export"

echo "RUNNING: qiime import DNASequence"
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path  $RAW/export_seqs/dna-sequences.fasta \
    --output-path $DB/ref-otus.qza
echo "DONE: qiime import DNASequence"

echo "RUNNING: qiime import Taxonomy"
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-format TSVTaxonomyFormat \
    --input-path  $RAW/export_tax/taxonomy.tsv \
    --output-path $DB/ref-taxonomy.qza
echo "DONE: qiime import Taxonomy"

echo "Outputs: $DB/ref-otus.qza and $DB/ref-taxonomy.qza"
}
