#!/usr/bin/env bash
# qiime2.sh - per-target QIIME2 amplicon pipeline (holoruminant).
#
# Author: Jeremy Tournayre
#
# Runs the full per-target workflow: primer check, import, primer/adapter
# trimming, DADA2 denoising, optional chimera removal, taxonomy, phylogeny,
# rarefaction, feature/sample filtering, then per-comparison diversity and
# differential abundance.
#
# Analysis parameters are read from a parameter file; all machine-specific
# paths come from config.sh.
#
# Usage (called by launch.pl):
#   bash qiime2.sh <TARGET_DIR> <RUN_CONF>
#
# The QIIME2 environment used for this run is CONDA_ENV_QIIME2 as set in
# config.sh. Before rerunning an existing project under a different QIIME2
# version, make sure config.sh reflects the version you intend to use.
#
# The pipeline is resumable: every step is guarded by `[[ ! -f output ]]`, and
# a few sub-steps intentionally continue after a failure (empty-.qzv fallback).
{

# --- Load central configuration and activate the QIIME2 environment ----------
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
activate_qiime2

# --- Inputs ------------------------------------------------------------------
# Resolved to an ABSOLUTE path immediately: TARGET_DIR is used to build
# TMPDIR below, which stays exported (and re-read by subprocesses) long
# after this script `cd`s into $TARGET_DIR/QIIME2. A relative TARGET_DIR
# would make TMPDIR resolve relative to the NEW working directory once
# subprocesses read it, doubling up the path (e.g.
# .../QIIME2/Data/<env>/<target>/tmp/... instead of .../<target>/tmp/...).
TARGET_DIR="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"   # Data/<target>, absolute
RUN_CONF="$2"       # per-target parameter file written by launch.pl

# --- Record pipeline provenance alongside this run's outputs -----------------
# Writes the exact git commit (and any uncommitted changes) of the pipeline
# code used for this run, so results can always be traced back to the code
# that produced them.
PROVENANCE_FILE="$TARGET_DIR/PROVENANCE.txt"
{
    echo "Pipeline run started: $(date -Iseconds)"
    echo "Pipeline directory:   $PIPELINE_DIR"
    echo "QIIME2 environment:   $CONDA_ENV_QIIME2 (from config.sh)"
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

# --- Load the named analysis parameters (KEY=value lines) --------------------
# shellcheck disable=SC1090
source "$RUN_CONF"
NSLOTS="${NSLOTS:-4}"
# Tag combining the database and its primer set, used to name taxonomy outputs.
TAX_PREFIX="${TAXONOMY_DATABASE}-${TAXONOMY_DATABASE_PRIMER}"

# --- Private per-run temp dir -----------------------------------------------
RANDOM_STRING=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 12)
export TMPDIR="$TARGET_DIR/tmp/$RANDOM_STRING"
mkdir -p "$TMPDIR"

cd "$TARGET_DIR/QIIME2" || exit 1

echo "------------------------------------------"
echo "target             $TARGET_DIR"
echo "run_conf           $RUN_CONF"
echo "pipeline_dir       $PIPELINE_DIR"
echo "primer_F           $PRIMER_F"
echo "primer_R           $PRIMER_R"
echo "primer_Rreverse    $PRIMER_RREVERSE"
echo "primer_Freverse    $PRIMER_FREVERSE"
echo "remove_reverse     $REMOVE_REVERSE"
echo "trim_F             $TRIM_F"
echo "trim_R             $TRIM_R"
echo "verif_primer_file  $VERIF_PRIMER_FILE"
echo "p-min-overlap      $P_MIN_OVERLAP"
echo "taxonomy_database  $TAXONOMY_DATABASE"
echo "taxonomy_db_primer $TAXONOMY_DATABASE_PRIMER"
echo "comparison         $COMPARISON"
echo "filter_name        $FILTER_NAME"
echo "filter             $FILTER"
echo "sampling_depth     $SAMPLING_DEPTH"
echo "single_or_paired   $SINGLE_OR_PAIRED"
echo "filtFeatures_mode  $FILTFEATURES_MODE"
echo "filtFeatures_val   $FILTFEATURES_VAL"
echo "boxplot            $BOXPLOT"
echo "trunc_F            $TRUNC_F"
echo "trunc_R            $TRUNC_R"
echo "rarefaction        $RAREFACTION"
echo "only_all_group     $ONLY_ALL_GROUP"
echo "p-min-fold-parent  $P_MIN_FOLD_PARENT"
echo "vsearch            $VSEARCH"
echo "ANCOMBC_Nb_Perc    $ANCOMBC_NB_PERC"

filter="$FILTER"

# --- Primer N-masking for the zgrep primer check ----------------------------
# reverse_f = primer_Rreverse, front_f = primer_F,
# reverse_r = primer_Freverse, front_r = primer_R.
secondString="."
reverse_f="$PRIMER_RREVERSE"; reverse_f_N="${reverse_f//[^ATCG]/"$secondString"}"
front_f="$PRIMER_F";          front_f_N="${front_f//[^ATCG]/"$secondString"}"
reverse_r="$PRIMER_FREVERSE"; reverse_r_N="${reverse_r//[^ATCG]/"$secondString"}"
front_r="$PRIMER_R";          front_r_N="${front_r//[^ATCG]/"$secondString"}"
seq_file_R1=`ls -1 ../input/*"$VERIF_PRIMER_FILE"*R1* | head -1`
seq_file_R2=`ls -1 ../input/*"$VERIF_PRIMER_FILE"*R2* | head -1`

verification_result () {
	if [[ ! -f $1 ]] && [[ ! -d  $1 ]]
	then
		echo "ERROR: the command above gave no result"
		exit
	fi
}

if [[ ! -f verif_primers_input.txt ]]
then
	echo "RUNNING: zgrep primers verification input file"
	echo "zgrep reverse_f_N seq_file_R1 $reverse_f_N $seq_file_R1" > verif_primers_input.txt
	zgrep $reverse_f_N $seq_file_R1 | wc -l  >> verif_primers_input.txt
	echo "zgrep front_f_N seq_file_R1 $front_f_N $seq_file_R1" >> verif_primers_input.txt
	zgrep $front_f_N $seq_file_R1 | wc -l >> verif_primers_input.txt
	echo "zgrep reverse_r_N seq_file_R1 $reverse_r_N $seq_file_R1" >> verif_primers_input.txt
	zgrep $reverse_r_N $seq_file_R1 | wc -l >> verif_primers_input.txt
	echo "zgrep front_r_N seq_file_R1 $front_r_N $seq_file_R1" >> verif_primers_input.txt
	zgrep $front_r_N $seq_file_R1 | wc -l >> verif_primers_input.txt

	echo "zgrep reverse_f_N seq_file_R2 $reverse_r_N $seq_file_R2" >> verif_primers_input.txt
	zgrep $reverse_f_N $seq_file_R2 | wc -l >> verif_primers_input.txt
	echo "zgrep front_f_N seq_file_R2 $front_r_N $seq_file_R2" >> verif_primers_input.txt
	zgrep $front_f_N $seq_file_R2 | wc -l >> verif_primers_input.txt
	echo "zgrep reverse_r_N seq_file_R2 $reverse_r_N $seq_file_R2" >> verif_primers_input.txt
	zgrep $reverse_r_N $seq_file_R2 | wc -l >> verif_primers_input.txt
	echo "zgrep front_r_N seq_file_R2 $front_r_N $seq_file_R2" >> verif_primers_input.txt
	zgrep $front_r_N $seq_file_R2 | wc -l >> verif_primers_input.txt
	echo "DONE: zgrep primers verification input file"
	verification_result "verif_primers_input.txt"
fi

if [[ ! -f QC/multiqc_report.html ]]
then
    echo "RUNNING: FastQC + MultiQC on raw reads"
    mkdir -p QC
    fastqc \
        --threads $NSLOTS \
        --outdir QC \
        input/*.fastq.gz
    multiqc \
        --outdir QC \
        --filename multiqc_report \
        --force \
        QC
    echo "DONE: FastQC + MultiQC on raw reads"
    verification_result "QC/multiqc_report.html"
fi

file="demux-paired-end"
mode="PairedEndSequencesWithQuality"
if [[ $SINGLE_OR_PAIRED == "single" ]]
then
	file="demux-single-end"
	mode="SequencesWithQuality"
fi

if [[ ! -f "$file.qza" ]] && [[ ! -f trimmed_remove_primers_wild.qza ]]
then
	echo "RUNNING: qiime import $mode"
	qiime tools import \
	  --type "SampleData[$mode]" \
	  --input-path ../input \
	  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
	  --output-path "$file.qza"
	echo "DONE: qiime import $mode"
	verification_result "$file.qza"
fi

if [[ ! -f "$file.qzv" ]]
then
	echo "RUNNING: qiime demux summarize"
	qiime demux summarize \
	  --i-data "$file.qza" \
	  --o-visualization "$file.qzv"
	echo "DONE: qiime demux summarize"
	verification_result "$file.qzv"
fi

input_file="$file.qza"
if [[ $PRIMER_F != "no" ]]
then
	if [ $REMOVE_REVERSE != "no" ]
	then
		if [[ ! -f trimmed_remove_adapter_wild.qza ]] &&  [[ ! -f trimmed_remove_primers_wild.qza ]]
		then
			if [ $SINGLE_OR_PAIRED == "single" ]
			then
				echo "RUNNING: qiime cutadapt trim-single trimmed_remove_adapter_wild.qza"
				qiime cutadapt trim-single \
				--i-demultiplexed-sequences "$file.qza" \
				--p-adapter $PRIMER_RREVERSE  \
				--p-cores $NSLOTS \
				--p-match-read-wildcards \
				--p-match-adapter-wildcards \
				--o-stats trimmed_remove_adapter_wild_stats.qza \
				--o-trimmed-sequences trimmed_remove_adapter_wild.qza 
				echo "DONE: qiime cutadapt trim-single trimmed_remove_adapter_wild.qza"
				verification_result "trimmed_remove_adapter_wild.qza"
				rm "$file.qza"
			else
				echo "RUNNING: qiime cutadapt trim-paired trimmed_remove_adapter_wild.qza"
				qiime cutadapt trim-paired \
				--i-demultiplexed-sequences "$file.qza" \
				--p-adapter-f $PRIMER_RREVERSE  \
				--p-adapter-r $PRIMER_FREVERSE \
				--p-cores $NSLOTS \
				--p-match-read-wildcards \
				--p-match-adapter-wildcards \
				--o-stats trimmed_remove_adapter_wild_stats.qza \
				--o-trimmed-sequences trimmed_remove_adapter_wild.qza 
				echo "DONE: qiime cutadapt trim-paired trimmed_remove_adapter_wild.qza"
				verification_result "trimmed_remove_adapter_wild.qza"
				rm "$file.qza"
			fi
		fi
		input_file="trimmed_remove_adapter_wild.qza"
	fi

	if [[ ! -f trimmed_remove_primers_wild.qza ]]
	then
		if [ $SINGLE_OR_PAIRED == "single" ]
		then
			echo "RUNNING: qiime cutadapt trim-single trimmed_remove_primers_wild.qza"
			qiime cutadapt trim-single \
			--i-demultiplexed-sequences $input_file \
			--p-front $PRIMER_F  \
			--p-match-read-wildcards \
			--p-cores $NSLOTS \
			--p-match-adapter-wildcards \
			--p-discard-untrimmed \
			--o-stats trimmed_remove_primers_wild_stats.qza \
			--o-trimmed-sequences trimmed_remove_primers_wild.qza 
			echo "DONE: qiime cutadapt trim-single trimmed_remove_primers_wild.qza"
			verification_result "trimmed_remove_primers_wild.qza"
			rm $input_file
		else
			echo "RUNNING: qiime cutadapt trim-paired trimmed_remove_primers_wild.qza"
			qiime cutadapt trim-paired \
			--i-demultiplexed-sequences $input_file \
			--p-front-f $PRIMER_F  \
			--p-front-r $PRIMER_R \
			--p-match-read-wildcards \
			--p-match-adapter-wildcards \
			--p-discard-untrimmed \
			--o-stats trimmed_remove_adapter_wild_stats.qza \
			--o-trimmed-sequences trimmed_remove_primers_wild.qza 
			echo "DONE: qiime cutadapt trim-paired trimmed_remove_primers_wild.qza"
			verification_result "trimmed_remove_primers_wild.qza"
			rm $input_file
		fi
	fi
	input_file="trimmed_remove_primers_wild.qza"

	if [[ ! -f verif_primers_trim.txt ]]
	then
		echo "RUNNING: qiime tools extract trimmed_remove_primers_wild"
		qiime tools extract --output-path trimmed_remove_primers_wild --input-path trimmed_remove_primers_wild.qza
		echo "DONE: qiime tools extract trimmed_remove_primers_wild"

		seq_dir=`ls trimmed_remove_primers_wild`
		seq_file_R1=`ls -1 trimmed_remove_primers_wild/$seq_dir/data/*"$VERIF_PRIMER_FILE"*R1* | head -1`
		seq_file_R2=`ls -1 trimmed_remove_primers_wild/$seq_dir/data/*"$VERIF_PRIMER_FILE"*R2* | head -1`

		echo "RUNNING: zgrep primers/reverses verification trimmed_remove_primers_wild"
		echo "zgrep reverse_f_N seq_file_R1 $reverse_f_N $seq_file_R1" > verif_primers_trim.txt
		zgrep $reverse_f_N $seq_file_R1 | wc -l  >> verif_primers_trim.txt
		echo "zgrep front_f_N seq_file_R1 $front_f_N $seq_file_R1" >> verif_primers_trim.txt
		zgrep $front_f_N $seq_file_R1 | wc -l >> verif_primers_trim.txt
		echo "zgrep reverse_r_N seq_file_R1 $reverse_r_N $seq_file_R1" >> verif_primers_trim.txt
		zgrep $reverse_r_N $seq_file_R1 | wc -l >> verif_primers_trim.txt
		echo "zgrep front_r_N seq_file_R1 $front_r_N $seq_file_R1" >> verif_primers_trim.txt
		zgrep $front_r_N $seq_file_R1 | wc -l >> verif_primers_trim.txt

		echo "zgrep reverse_f_N seq_file_R2 $reverse_r_N $seq_file_R2" >> verif_primers_trim.txt
		zgrep $reverse_f_N $seq_file_R2 | wc -l >> verif_primers_trim.txt
		echo "zgrep front_f_N seq_file_R2 $front_r_N $seq_file_R2" >> verif_primers_trim.txt
		zgrep $front_f_N $seq_file_R2 | wc -l >> verif_primers_trim.txt
		echo "zgrep reverse_r_N seq_file_R2 $reverse_r_N $seq_file_R2" >> verif_primers_trim.txt
		zgrep $reverse_r_N $seq_file_R2 | wc -l >> verif_primers_trim.txt
		echo "zgrep front_r_N seq_file_R2 $front_r_N $seq_file_R2" >> verif_primers_trim.txt
		zgrep $front_r_N $seq_file_R2 | wc -l >> verif_primers_trim.txt

		echo "DONE: zgrep primers/reverses verification trimmed_remove_primers_wild"
		rm -R trimmed_remove_primers_wild
		verification_result "verif_primers_trim.txt"
	fi

	if [[ ! -f trimmed_remove_primers_wild.qzv ]]
	then
		echo "RUNNING: qiime demux summarize"
		qiime demux summarize \
		--i-data trimmed_remove_primers_wild.qza \
		--o-visualization trimmed_remove_primers_wild.qzv
		echo "DONE: qiime demux summarize"
		verification_result "trimmed_remove_primers_wild.qzv"
	fi
fi

############## Stop here until trunc/trim parameters are defined from the .qzv
if [ "$TRIM_F" == "NA" ] || [ "$TRIM_R" == "NA" ] || [ "$TRUNC_F" == "NA" ] || [ "$TRUNC_R" == "NA" ]
then
	echo "EXIT: You can define trunc/trim parameter"
	exit
fi

if [[ ! -f rep-seqs.qza ]]
then
	if [ $SINGLE_OR_PAIRED == "single" ]
	then
		echo "RUNNING: qiime dada2 denoise-single"
		qiime dada2 denoise-single \
		--i-demultiplexed-seqs $input_file \
		--p-trim-left $TRIM_F \
		--p-trunc-len $TRUNC_F \
		--p-min-fold-parent-over-abundance $P_MIN_FOLD_PARENT \
		--o-table table.qza \
		--o-representative-sequences rep-seqs.qza \
		--o-base-transition-stats base-transition-stats.qza \
		--p-n-threads $NSLOTS \
		--o-denoising-stats denoising-stats.qza
		echo "DONE: qiime dada2 denoise-single"
		verification_result "rep-seqs.qza"
	else
		echo "RUNNING: qiime dada2 denoise-paired"
		qiime dada2 denoise-paired \
		--i-demultiplexed-seqs $input_file \
		--p-trim-left-f $TRIM_F \
		--p-trim-left-r $TRIM_R \
		--p-trunc-len-f $TRUNC_F \
		--p-trunc-len-r $TRUNC_R \
		--p-n-threads $NSLOTS \
		--p-min-fold-parent-over-abundance $P_MIN_FOLD_PARENT \
		--o-table table.qza \
		--p-min-overlap $P_MIN_OVERLAP \
		--o-representative-sequences rep-seqs.qza \
		--o-base-transition-stats base-transition-stats.qza \
		--o-denoising-stats denoising-stats.qza
		echo "DONE: qiime dada2 denoise-paired"	
		verification_result "rep-seqs.qza"
	fi
fi

if [[ ! -f denoising-stats.qzv ]]
then
	echo "RUNNING: visualization denoising-stats.qzv"
	qiime metadata tabulate \
	--m-input-file denoising-stats.qza \
	--o-visualization denoising-stats.qzv
	echo "DONE: visualization denoising-stats.qzv"
	verification_result "denoising-stats.qzv"
fi

if [[ ! -f ../mapping.tsv ]]
then
	echo "EXIT: You must provide mapping.tsv file, you can add depth data at this step"
	exit 
fi

if [[ ! -f table.qzv ]]
then
	echo "RUNNING: visualization table.qzv"
	qiime feature-table summarize \
	--i-table table.qza \
	--m-metadata-file ../mapping.tsv \
	--o-feature-frequencies feature-frequencies.qza \
	--o-sample-frequencies sample-frequencies.qza \
	--o-summary table.qzv
	echo "DONE: visualization table.qzv"
	verification_result "table.qzv"
fi

if [[ ! -f rep-seqs.qzv ]]
then
	echo "RUNNING: visualization rep-seqs.qzv"
	qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv
	echo "DONE: visualization rep-seqs.qzv"
	verification_result "rep-seqs.qzv"
fi


dir_table="./"
dir_table_rep_seqs="./"
if [[ $VSEARCH == "vsearch" ]]
then
	dir_table="./vsearch_uchime-denovo/"
	dir_table_rep_seqs="./vsearch_uchime-denovo/"
	if [[ ! -f "vsearch_uchime-denovo/stats.qza" ]]
	then
		echo "RUN: vsearch_uchime-denovo"
		qiime vsearch uchime-denovo \
		--i-sequences rep-seqs.qza \
		--i-table table.qza \
		--p-mindiv 0.2 \
		--p-minh 0.05 \
		--output-dir vsearch_uchime-denovo 
		echo "DONE: vsearch_uchime-denovo"
		verification_result "vsearch_uchime-denovo/stats.qza"
	fi

	if [[ ! -f "vsearch_uchime-denovo/stats.qzv" ]]
	then
		echo "RUN: vsearch visualization stats.qzv"
		qiime metadata tabulate \
		--m-input-file vsearch_uchime-denovo/stats.qza \
		--o-visualization vsearch_uchime-denovo/stats.qzv
		echo "DONE: vsearch visualization stats.qzv"
		verification_result "vsearch_uchime-denovo/stats.qzv"
	fi

	if [[ ! -f "vsearch_uchime-denovo/table.qza" ]]
	then
		echo "RUN: vsearch filter-features table.qza"
		qiime feature-table filter-features \
		--i-table table.qza \
		--m-metadata-file vsearch_uchime-denovo/nonchimeras.qza \
		--o-filtered-table vsearch_uchime-denovo/table.qza
		echo "DONE: vsearch filter-features table.qza"
		verification_result "vsearch_uchime-denovo/table.qza"
	fi

	if [[ ! -f "vsearch_uchime-denovo/rep-seqs.qza" ]]
	then
		echo "RUN: vsearch filter-seqs rep-seqs.qza"
		qiime feature-table filter-seqs \
		--i-data rep-seqs.qza \
		--m-metadata-file vsearch_uchime-denovo/nonchimeras.qza \
		--o-filtered-data vsearch_uchime-denovo/rep-seqs.qza
		echo "DONE: vsearch filter-seqs rep-seqs.qza"
		verification_result "vsearch_uchime-denovo/rep-seqs.qza"
	fi

	if [[ ! -f "vsearch_uchime-denovo/table.qzv" ]]
	then
		echo "RUNNING: visualization vsearch_uchime-denovo/table.qzv"
		qiime feature-table summarize \
		--i-table vsearch_uchime-denovo/table.qza \
		--m-metadata-file ../mapping.tsv \
		--o-feature-frequencies vsearch_uchime-denovo/feature-frequencies.qza \
		--o-sample-frequencies vsearch_uchime-denovo/sample-frequencies.qza \
		--o-summary vsearch_uchime-denovo/table.qzv
		echo "DONE: visualization vsearch_uchime-denovo/table.qzv"
		verification_result "vsearch_uchime-denovo/table.qzv"
	fi

	if [[ ! -f "vsearch_uchime-denovo/rep-seqs.qzv" ]]
	then
		echo "RUNNING: visualization vsearch_uchime-denovo/rep-seqs.qzv"
		qiime feature-table tabulate-seqs \
		--i-data vsearch_uchime-denovo/rep-seqs.qza \
		--o-visualization vsearch_uchime-denovo/rep-seqs.qzv
		echo "DONE: visualization vsearch_uchime-denovo/rep-seqs.qzv"
		verification_result "vsearch_uchime-denovo/rep-seqs.qzv"
	fi
fi


if [[ ! -f $dir_table/$TAX_PREFIX-taxonomy.qza ]]
then
	echo "RUNNING: classification"
	qiime feature-classifier classify-sklearn \
	--i-classifier $DB_DIR/Data/$TAXONOMY_DATABASE/$CONDA_ENV_QIIME2/$TAXONOMY_DATABASE_PRIMER/classifier.qza \
	--i-reads $dir_table/rep-seqs.qza \
	--p-n-jobs $NSLOTS \
	--o-classification $dir_table/$TAX_PREFIX-taxonomy.qza
	echo "DONE: classification"
	verification_result "$dir_table/$TAX_PREFIX-taxonomy.qza"
fi

if [[ ! -f $dir_table/$TAX_PREFIX-taxonomy.qzv ]]
then
	echo "RUNNING: taxonomy"
	qiime metadata tabulate \
	--m-input-file $dir_table/$TAX_PREFIX-taxonomy.qza \
	--o-visualization $dir_table/$TAX_PREFIX-taxonomy.qzv
	echo "DONE: taxonomy"
	verification_result "$dir_table/$TAX_PREFIX-taxonomy.qzv"
fi


if [[ ! -f $dir_table/$TAX_PREFIX-taxa-bar-plots.qzv ]]
then
	echo "RUNNING:  taxa barplot"  
	qiime taxa barplot \
	--i-table $dir_table/table.qza \
	--i-taxonomy $dir_table/$TAX_PREFIX-taxonomy.qza \
	--m-metadata-file ../mapping.tsv \
	--o-visualization $dir_table/$TAX_PREFIX-taxa-bar-plots.qzv
	echo "DONE:  taxa barplot"
	verification_result "$dir_table/$TAX_PREFIX-taxa-bar-plots.qzv"
fi

if [[ ! -f $dir_table/$TAX_PREFIX-taxa-bar-plots.xlsx ]]
then
	conda deactivate
	echo "Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_table/table.qza $dir_table/$TAX_PREFIX-taxonomy.qza $dir_table/$TAX_PREFIX-taxa-bar-plots.xlsx"
	echo "RUNNING: Rscript extract_taxo_barplots"
	Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_table/table.qza $dir_table/$TAX_PREFIX-taxonomy.qza $dir_table/$TAX_PREFIX-taxa-bar-plots.xlsx
	echo "DONE: Rscript extract_taxo_barplots"
	activate_qiime2
fi

if [[ ! -f $dir_table/rooted-tree.qza ]]
then
	echo "RUNNING: qiime phylogeny"
	qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences $dir_table/rep-seqs.qza \
	--o-alignment $dir_table/aligned-rep-seqs.qza \
	--p-n-threads $NSLOTS \
	--o-masked-alignment $dir_table/masked-aligned-rep-seqs.qza \
	--o-tree $dir_table/unrooted-tree.qza \
	--o-rooted-tree $dir_table/rooted-tree.qza
	echo "DONE: qiime phylogeny"
	verification_result "$dir_table/rooted-tree.qza"
fi
if ! ls $dir_table/alpha-rarefaction_*.qzv 1>/dev/null 2>&1
then
	if [[ ! -d $dir_table/table ]]
	then
		cp $dir_table/table.qza $dir_table/table.zip
		unzip -q $dir_table/table.zip -d $dir_table/table
	fi
	code=`ls  $dir_table/table`;
	conda deactivate

	`Rscript $RSCRIPTS_DIR/extract_samples.R $dir_table $code`

	rarefaction_depth=`Rscript $RSCRIPTS_DIR/get_lowest_depth.R $dir_table $code`
	rarefaction_depth_max=`Rscript $RSCRIPTS_DIR/get_max_depth.R $dir_table $code`
	rarefaction_depth="${rarefaction_depth:4}"
	rarefaction_depth_max="${rarefaction_depth_max:4}"
	curve_rarefaction=$rarefaction_depth
	activate_qiime2
	if (( rarefaction_depth_max < curve_rarefaction ))
	then
		curve_rarefaction=$rarefaction_depth_max
	fi

	echo "RUNNING: alpha-rarefaction $dir_table"
	qiime diversity alpha-rarefaction \
	--i-table $dir_table/table.qza \
	--i-phylogeny  $dir_table/rooted-tree.qza \
	--p-max-depth $curve_rarefaction \
	--p-steps 5 \
	--p-iterations 1 \
	--m-metadata-file  ../mapping.tsv \
	--o-visualization $dir_table/alpha-rarefaction_$curve_rarefaction.qzv
	echo "DONE: alpha-rarefaction $dir_table"
	verification_result "$dir_table/alpha-rarefaction_$curve_rarefaction.qzv"
fi

if [[ ! -f $dir_table/goodcov_Mean.txt ]]
then
	echo "RUNNING: Rscript goodcoverage_report"
	conda deactivate
	echo "Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_table $dir_table/table.qza"
	Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_table $dir_table/table.qza
	echo "DONE: Rscript goodcoverage_report"
	activate_qiime2
fi	


dir_filt_feature="Filter_Stat/All"
##########################################################
if [[ $FILTFEATURES_MODE != "All" ]]
then
	dir_filt_feature="Filter_Stat/$FILTFEATURES_MODE-$FILTFEATURES_VAL"
	if [[ $FILTFEATURES_MODE == "minSample" ]] && [[ ! -f $dir_filt_feature/table.qza ]]
	then
		if [[ ! -d  $dir_filt_feature ]]
		then
			mkdir -p $dir_filt_feature
		fi		
		echo "RUNNING:  feature-table filter-features"
		qiime feature-table filter-features \
		  --i-table $dir_table/table.qza \
		  --p-min-samples $FILTFEATURES_VAL \
		  --o-filtered-table $dir_filt_feature/table.qza
		 echo "DONE:  feature-table filter-features"
		 verification_result "$dir_filt_feature/table.qza"
	elif [[ $FILTFEATURES_MODE == "file" ]]
	then
		dir_filt_feature="Filter_Stat/$FILTFEATURES_MODE-$FILTFEATURES_VAL-$TAX_PREFIX"
		if [[ ! -f $dir_filt_feature/table.qza ]]
		then
			if [[ ! -d  $dir_filt_feature ]]
			then
				mkdir -p $dir_filt_feature
			fi	
			echo "RUNNING:  taxa filter-table"
			filt_taxo=`cat ../Filter/Features/$FILTFEATURES_VAL`
			qiime taxa filter-table \
			  --i-table $dir_table/table.qza \
			  --i-taxonomy $dir_table/$TAX_PREFIX-taxonomy.qza \
			  --p-include $filt_taxo \
			  --o-filtered-table $dir_filt_feature/table.qza
			echo "DONE:  taxa filter-table"	
			verification_result "$dir_filt_feature/table.qza"
		fi
	fi
	if [[ ! -f $dir_filt_feature/transposed_table.qza ]]
	then
		echo "RUNNING:  feature-table transpose"
		qiime feature-table transpose \
		 --i-table $dir_filt_feature/table.qza \
		 --o-transposed-feature-table $dir_filt_feature/transposed_table.qza
		echo "DONE:  feature-table transpose"
		verification_result "$dir_filt_feature/transposed_table.qza"
	fi

	if [[ ! -f $dir_filt_feature/$TAX_PREFIX-taxonomy.qza ]]
	then
		echo "RUNNING:  rescript filter-taxa"  
		qiime rescript filter-taxa \
		--i-taxonomy $dir_table/$TAX_PREFIX-taxonomy.qza \
		--m-ids-to-keep-file $dir_filt_feature/transposed_table.qza \
		--o-filtered-taxonomy $dir_filt_feature/$TAX_PREFIX-taxonomy.qza
		echo "DONE:  rescript filter-taxa"
		verification_result "$dir_filt_feature/$TAX_PREFIX-taxonomy.qza"
	fi

	if [[ ! -f $dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.qzv ]]
	then
		echo "RUNNING:  taxa barplot"  
		qiime taxa barplot \
		--i-table $dir_filt_feature/table.qza \
		--i-taxonomy $dir_filt_feature/$TAX_PREFIX-taxonomy.qza \
		--m-metadata-file ../mapping.tsv \
		--o-visualization $dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.qzv
		echo "DONE:  taxa barplot"
		verification_result "$dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.qzv"
	fi
	if [[ ! -f $dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.xlsx ]]
	then
		conda deactivate
		echo "Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_filt_feature/table.qza $dir_filt_feature/$TAX_PREFIX-taxonomy.qza $dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.xlsx"
		echo "RUNNING: Rscript extract_taxo_barplots"
		Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_filt_feature/table.qza $dir_filt_feature/$TAX_PREFIX-taxonomy.qza $dir_filt_feature/$TAX_PREFIX-taxa-bar-plots.xlsx
		echo "DONE: Rscript extract_taxo_barplots"
		activate_qiime2
	fi	
	dir_table=$dir_filt_feature
else
	if [[ ! -d  $dir_filt_feature ]]
	then
		mkdir $dir_filt_feature
	fi	
fi


dir_filt_feature="$dir_filt_feature/$FILTER_NAME"
if [[ ! -d  $dir_filt_feature ]]
then
	mkdir $dir_filt_feature
fi	

if [[ ! -f ../Filter/Samples/$FILTER_NAME.txt ]]
then
	echo "RUNNING: perl filter mapping $FILTER_NAME $filter"   
	echo "perl $PLSCRIPTS_DIR/filter_mapping.pl $FILTER_NAME $filter"
	eval "perl $PLSCRIPTS_DIR/filter_mapping.pl $FILTER_NAME \"$filter\""
	echo "DONE: perl filter mapping $FILTER_NAME"
	verification_result "../Filter/Samples/$FILTER_NAME.txt"
fi

declare -A sampling_depth
IFS=','
read -a strarr <<< "$SAMPLING_DEPTH"
IFS=':'
for ANCOM_vs in "${strarr[@]}"
do
	read -a strarr2 <<< "$ANCOM_vs"
	sampling_depth[${strarr2[0]}]=${strarr2[1]}
done

unset IFS
dir_filt_feature="$dir_filt_feature/$COMPARISON"
if [[ ! -d  $dir_filt_feature ]]
then
	mkdir $dir_filt_feature
fi


#Pairwise
echo "perl $PLSCRIPTS_DIR/get_value_ancom.pl $FILTER_NAME $COMPARISON $dir_filt_feature $ONLY_ALL_GROUP"
pairwise=`perl $PLSCRIPTS_DIR/get_value_ancom.pl $FILTER_NAME $COMPARISON $dir_filt_feature $ONLY_ALL_GROUP`

#add ALL
if [[ ! -d  $dir_filt_feature/All ]]
then
	mkdir $dir_filt_feature/All
fi
cp ../Filter/Samples/$FILTER_NAME.txt $dir_filt_feature/All/samples.txt

save_dir_filt_feature=$dir_filt_feature
for ANCOM_vs in $pairwise
do
	dir_filt_feature="$save_dir_filt_feature/$ANCOM_vs"
	if [[ ! -f $dir_filt_feature/table.qza ]]
	then
		echo "RUNNING: feature-table filter-samples $dir_filt_feature"
		qiime feature-table filter-samples \
		--i-table $dir_table/table.qza \
		--m-metadata-file  $dir_filt_feature/samples.txt \
		--o-filtered-table $dir_filt_feature/table.qza
		echo "DONE: feature-table filter-samples $dir_filt_feature"
		verification_result "$dir_filt_feature/table.qza"
	fi 

	if [[ ! -f $dir_filt_feature/table.qzv ]]
	then
		echo "RUNNING: feature-table summarize $dir_filt_feature"
		qiime feature-table summarize \
		--i-table $dir_filt_feature/table.qza \
		--m-metadata-file ../mapping.tsv \
		--o-feature-frequencies $dir_filt_feature/feature-frequencies.qza \
		--o-sample-frequencies $dir_filt_feature/sample-frequencies.qza \
		--o-summary $dir_filt_feature/table.qzv
		echo "DONE: feature-table summarize $dir_filt_feature"
		verification_result "$dir_filt_feature/table.qzv"
	fi	

	
	if [[ $RAREFACTION == "yes" ]]
	then 
		if [[ ! -f  $dir_filt_feature/rep-seqs.qza ]]
		then		
			echo "RUNNING: feature-table filter-seqs $dir_filt_feature"
			qiime feature-table filter-seqs \
			--i-data $dir_table_rep_seqs/rep-seqs.qza \
			--i-table $dir_filt_feature/table.qza \
			--o-filtered-data  $dir_filt_feature/rep-seqs.qza
			echo "DONE: feature-table filter-seqs $dir_filt_feature"
			verification_result "$dir_filt_feature/rep-seqs.qza"
		fi	
		if [[ ! -f $dir_filt_feature/rooted-tree.qza ]]
		then
			echo "RUNNING: qiime phylogeny $dir_filt_feature"
			qiime phylogeny align-to-tree-mafft-fasttree \
			--i-sequences $dir_filt_feature/rep-seqs.qza \
			--o-alignment $dir_filt_feature/aligned-rep-seqs.qza \
			--p-n-threads $NSLOTS \
			--o-masked-alignment $dir_filt_feature/masked-aligned-rep-seqs.qza \
			--o-tree $dir_filt_feature/unrooted-tree.qza \
			--o-rooted-tree $dir_filt_feature/rooted-tree.qza
			echo "DONE: qiime phylogeny $dir_filt_feature"
			verification_result "$dir_filt_feature/rooted-tree.qza"
		fi
		if ! ls $dir_filt_feature/alpha-rarefaction_*.qzv 1>/dev/null 2>&1
		then
			if [[ ! -d $dir_filt_feature/table ]]
			then
				cp $dir_filt_feature/table.qza $dir_filt_feature/table.zip
				unzip -q $dir_filt_feature/table.zip -d $dir_filt_feature/table
			fi
			code=`ls  $dir_filt_feature/table`;
			conda deactivate

			`Rscript $RSCRIPTS_DIR/extract_samples.R $dir_filt_feature $code`

			rarefaction_depth=`Rscript $RSCRIPTS_DIR/get_lowest_depth.R $dir_filt_feature $code`
			rarefaction_depth_max=`Rscript $RSCRIPTS_DIR/get_max_depth.R $dir_filt_feature $code`
			rarefaction_depth="${rarefaction_depth:4}"
			rarefaction_depth_max="${rarefaction_depth_max:4}"
			curve_rarefaction=$rarefaction_depth
			activate_qiime2
			if (( rarefaction_depth_max < curve_rarefaction ))
			then
				curve_rarefaction=$rarefaction_depth_max
			fi

			echo "RUNNING: alpha-rarefaction $dir_filt_feature"
			qiime diversity alpha-rarefaction \
			--i-table $dir_filt_feature/table.qza \
			--i-phylogeny  $dir_filt_feature/rooted-tree.qza \
			--p-max-depth $curve_rarefaction \
			--p-steps 5 \
			--p-iterations 1 \
			--m-metadata-file  ../mapping.tsv \
			--o-visualization $dir_filt_feature/alpha-rarefaction_$curve_rarefaction.qzv
			echo "DONE: alpha-rarefaction $dir_filt_feature"
			verification_result "$dir_filt_feature/alpha-rarefaction_$curve_rarefaction.qzv"
		fi

		if [[ ! -f $dir_filt_feature/goodcov_Mean.txt ]]
		then
			echo "RUNNING: Rscript goodcoverage_report"
			conda deactivate
			echo "Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_filt_feature $dir_filt_feature/table.qza"
			Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_filt_feature $dir_filt_feature/table.qza
			echo "DONE: Rscript goodcoverage_report"
			activate_qiime2
		fi		
	fi
	if [[ ${sampling_depth[$ANCOM_vs]} ]]
	then
		dir_filt_feature=$dir_filt_feature/${sampling_depth[$ANCOM_vs]}
		if [[ ! -d  $dir_filt_feature ]]
		then
			mkdir $dir_filt_feature
		fi
		if [[ ! -f  $dir_filt_feature/table.qza ]] &&  [[ ${sampling_depth[$ANCOM_vs]} != 0 ]]
		then				
			echo "RUNNING: feature-table filter-samples $dir_filt_feature/table.qza"
			qiime feature-table filter-samples \
			--i-table $dir_filt_feature/../table.qza \
			--p-min-frequency ${sampling_depth[$ANCOM_vs]} \
			--o-filtered-table  $dir_filt_feature/table.qza
			echo "DONE: feature-table filter-samples $dir_filt_feature/table.qza"
			verification_result "$dir_filt_feature/table.qza"
		elif [[ ${sampling_depth[$ANCOM_vs]} == 0 ]]
		then
			cp $dir_filt_feature/../table.qza $dir_filt_feature/table.qza
		fi

		if [[ ! -f  $dir_filt_feature/rep-seqs.qza ]]
		then		
			echo "RUNNING: feature-table filter-seqs $ANCOM_vs"
			qiime feature-table filter-seqs \
			--i-data $dir_table_rep_seqs/rep-seqs.qza \
			--i-table $dir_filt_feature/table.qza \
			--o-filtered-data  $dir_filt_feature/rep-seqs.qza
			echo "DONE: feature-table filter-seqs $ANCOM_vs"
			verification_result "$dir_filt_feature/rep-seqs.qza"
		fi	

		if [[ ! -f $dir_filt_feature/rooted-tree.qza ]]
		then
			echo "RUNNING: qiime phylogeny $dir_filt_feature"
			qiime phylogeny align-to-tree-mafft-fasttree \
			--i-sequences $dir_filt_feature/rep-seqs.qza \
			--o-alignment $dir_filt_feature/aligned-rep-seqs.qza \
			--p-n-threads $NSLOTS \
			--o-masked-alignment $dir_filt_feature/masked-aligned-rep-seqs.qza \
			--o-tree $dir_filt_feature/unrooted-tree.qza \
			--o-rooted-tree $dir_filt_feature/rooted-tree.qza
			echo "DONE: qiime phylogeny $dir_filt_feature"
			verification_result "$dir_filt_feature/rooted-tree.qza"
		fi

		########rarefaction to the lowest depth
		if [[ ! -d $dir_filt_feature/table ]]
		then
			cp $dir_filt_feature/table.qza $dir_filt_feature/table.zip
			unzip -q $dir_filt_feature/table.zip -d $dir_filt_feature/table
		fi
		code=`ls  $dir_filt_feature/table`;
		conda deactivate

		`Rscript $RSCRIPTS_DIR/extract_samples.R $dir_filt_feature $code`

		rarefaction_depth=`Rscript $RSCRIPTS_DIR/get_lowest_depth.R $dir_filt_feature $code`
		rarefaction_depth_max=`Rscript $RSCRIPTS_DIR/get_max_depth.R $dir_filt_feature $code`
		rarefaction_depth="${rarefaction_depth:4}"
		rarefaction_depth_max="${rarefaction_depth_max:4}"
		curve_rarefaction=$rarefaction_depth
		activate_qiime2
		if (( rarefaction_depth_max < curve_rarefaction ))
		then
			curve_rarefaction=$rarefaction_depth_max
		fi
		
		if ! ls $dir_filt_feature/alpha-rarefaction_*.qzv 1>/dev/null 2>&1
		then
			echo "RUNNING: alpha-rarefaction $dir_filt_feature"
			qiime diversity alpha-rarefaction \
			--i-table $dir_filt_feature/table.qza \
			--i-phylogeny  $dir_filt_feature/rooted-tree.qza \
			--p-max-depth $curve_rarefaction \
			--p-steps 5 \
			--p-iterations 1 \
			--m-metadata-file  ../mapping.tsv \
			--o-visualization $dir_filt_feature/alpha-rarefaction_$curve_rarefaction.qzv
			echo "DONE: alpha-rarefaction"
			verification_result "$dir_filt_feature/alpha-rarefaction_$curve_rarefaction.qzv"
		fi

		if [[ ! -f $dir_filt_feature/goodcov_Mean.txt ]]
		then
			echo "RUNNING: Rscript goodcoverage_report"
			conda deactivate
			echo "Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_filt_feature $dir_filt_feature/table.qza"
			Rscript $RSCRIPTS_DIR/goodcoverage_report.R $dir_filt_feature $dir_filt_feature/table.qza
			echo "DONE: Rscript goodcoverage_report"
			activate_qiime2
		fi		


		if [[ ! -d $dir_filt_feature/core-metrics-results_$rarefaction_depth ]]
		then

			echo "RUNNING: core-metrics-phylogenetic $dir_filt_feature $rarefaction_depth"
			qiime diversity core-metrics-phylogenetic \
				--i-phylogeny $dir_filt_feature/rooted-tree.qza \
				--i-table $dir_filt_feature/table.qza \
				--p-sampling-depth $rarefaction_depth \
				--p-n-jobs-or-threads $NSLOTS \
				--p-random-seed 42 \
				--m-metadata-file ../mapping.tsv \
				--output-dir $dir_filt_feature/core-metrics-results_$rarefaction_depth
			echo "DONE: core-metrics-phylogenetic $dir_filt_feature $rarefaction_depth"
			verification_result "$dir_filt_feature/core-metrics-results_$rarefaction_depth"
		fi
		save_dir_filt_feature_2=$dir_filt_feature
		dir_filt_feature=$dir_filt_feature/core-metrics-results_$rarefaction_depth

		for metric in weighted_unifrac unweighted_unifrac jaccard bray_curtis
		do
			if [[ ! -f $dir_filt_feature/metadata_${metric}.png ]]
			then
				echo "RUNNING: diversity ${metric}_pcoa_results Rplot $dir_filt_feature $rarefaction_depth"
				
				qiime tools export \
					--input-path $dir_filt_feature/${metric}_pcoa_results.qza \
					--output-path $dir_filt_feature/exported_pcoa_${metric}
				
				echo "DONE: diversity ${metric}_pcoa_results Rplot $dir_filt_feature $rarefaction_depth"

				conda deactivate

				echo "RUNNING: Rscript add_metadata_on_diversity_plot ${metric}_pcoa_results"
				echo "Rscript $RSCRIPTS_DIR/add_metadata_on_diversity_plot.R $dir_filt_feature/exported_pcoa_${metric}/ordination.txt ../mapping.tsv $COMPARISON $dir_filt_feature/metadata_${metric} ../mapping_extra.tsv"

				Rscript $RSCRIPTS_DIR/add_metadata_on_diversity_plot.R \
					$dir_filt_feature/exported_pcoa_${metric}/ordination.txt \
					../mapping.tsv \
					$COMPARISON \
					$dir_filt_feature/metadata_${metric} \
					../mapping_extra.tsv

				echo "DONE: Rscript add_metadata_on_diversity_plot ${metric}_pcoa_results"

				activate_qiime2
			fi
		done


		if [[ ! -f $dir_filt_feature/faith-pd-group-significance.qzv ]]
		then
			echo "RUNNING: diversity alpha-group-significance faith-pd $dir_filt_feature $rarefaction_depth"
			qiime diversity alpha-group-significance \
				--i-alpha-diversity $dir_filt_feature/faith_pd_vector.qza \
				--m-metadata-file ../mapping.tsv \
				--o-visualization $dir_filt_feature/faith-pd-group-significance.qzv
			if [[ ! -f $dir_filt_feature/faith-pd-group-significance.qzv ]]
			then
				touch $dir_filt_feature/faith-pd-group-significance.qzv
			fi				  
			echo "DONE: diversity alpha-group-significance faith-pd $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/evenness_vector-group-significance.qzv ]]
		then
			echo "RUNNING: diversity alpha-group-significance evenness_vector $dir_filt_feature $rarefaction_depth"
			qiime diversity alpha-group-significance \
				--i-alpha-diversity $dir_filt_feature/evenness_vector.qza \
				--m-metadata-file ../mapping.tsv \
				--o-visualization $dir_filt_feature/evenness_vector-group-significance.qzv
			if [[ ! -f $dir_filt_feature/evenness_vector-group-significance.qzv ]]
			then
				touch $dir_filt_feature/evenness_vector-group-significance.qzv
			fi					  
			echo "DONE: diversity alpha-group-significance evenness_vector $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/shannon_vector-group-significance.qzv ]]
		then
			echo "RUNNING: diversity alpha-group-significance shannon_vector $dir_filt_feature $rarefaction_depth"
			qiime diversity alpha-group-significance \
				--i-alpha-diversity $dir_filt_feature/shannon_vector.qza \
				--m-metadata-file ../mapping.tsv \
				--o-visualization $dir_filt_feature/shannon_vector-group-significance.qzv
			if [[ ! -f $dir_filt_feature/shannon_vector-group-significance.qzv ]]
			then
				touch $dir_filt_feature/shannon_vector-group-significance.qzv
			fi						  
			echo "DONE: diversity alpha-group-significance shannon_vector $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/observed_features_vector-group-significance.qzv ]]
		then
			echo "RUNNING: diversity alpha-group-significance observed_features_vector $dir_filt_feature $rarefaction_depth"
			qiime diversity alpha-group-significance \
				--i-alpha-diversity $dir_filt_feature/observed_features_vector.qza \
				--m-metadata-file ../mapping.tsv \
				--o-visualization $dir_filt_feature/observed_features_vector-group-significance.qzv
			if [[ ! -f $dir_filt_feature/observed_features_vector-group-significance.qzv ]]
			then
				touch $dir_filt_feature/observed_features_vector-group-significance.qzv
			fi						  
			echo "DONE: diversity alpha-group-significance observed_features_vector $dir_filt_feature $rarefaction_depth"
		fi

			if [[ ! -f $dir_filt_feature/weighted_unifrac_distance_matrix.qzv ]]
		then
			echo "RUNNING: diversity beta-group-significance weighted_unifrac_distance_matrix $dir_filt_feature $rarefaction_depth"
			qiime diversity beta-group-significance \
				--i-distance-matrix  $dir_filt_feature/weighted_unifrac_distance_matrix.qza \
				--m-metadata-file ../mapping.tsv \
				--m-metadata-column $COMPARISON \
				--o-visualization $dir_filt_feature/weighted_unifrac_distance_matrix.qzv \
				--p-pairwise
			if [[ ! -f $dir_filt_feature/weighted_unifrac_distance_matrix.qzv ]]
			then
				touch $dir_filt_feature/weighted_unifrac_distance_matrix.qzv
			fi						  
			echo "DONE: diversity beta-group-significance weighted_unifrac_distance_matrix $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/unweighted_unifrac_distance_matrix.qzv ]]
		then
			echo "RUNNING: diversity beta-group-significance unweighted_unifrac_distance_matrix $dir_filt_feature $rarefaction_depth"
			qiime diversity beta-group-significance \
				--i-distance-matrix  $dir_filt_feature/unweighted_unifrac_distance_matrix.qza \
				--m-metadata-file ../mapping.tsv \
				--m-metadata-column $COMPARISON \
				--o-visualization  $dir_filt_feature/unweighted_unifrac_distance_matrix.qzv \
				--p-pairwise
			if [[ ! -f $dir_filt_feature/unweighted_unifrac_distance_matrix.qzv ]]
			then
				touch $dir_filt_feature/unweighted_unifrac_distance_matrix.qzv
			fi						  
			echo "DONE: diversity beta-group-significance unweighted_unifrac_distance_matrix $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/jaccard_distance_matrix.qzv ]]
		then
			echo "RUNNING: diversity beta-group-significance jaccard_distance_matrix $dir_filt_feature $rarefaction_depth"
			qiime diversity beta-group-significance \
				--i-distance-matrix $dir_filt_feature/jaccard_distance_matrix.qza \
				--m-metadata-file ../mapping.tsv \
				--m-metadata-column $COMPARISON \
				--o-visualization $dir_filt_feature/jaccard_distance_matrix.qzv \
				--p-pairwise
			if [[ ! -f $dir_filt_feature/jaccard_distance_matrix.qzv ]]
			then
				touch $dir_filt_feature/jaccard_distance_matrix.qzv
			fi						  
			echo "DONE: diversity beta-group-significance jaccard_distance_matrix $dir_filt_feature $rarefaction_depth"
		fi

		if [[ ! -f $dir_filt_feature/bray_curtis_distance_matrix.qzv ]]
		then
			echo "RUNNING: diversity beta-group-significance bray_curtis_distance_matrix $dir_filt_feature $rarefaction_depth"
			qiime diversity beta-group-significance \
				--i-distance-matrix $dir_filt_feature/bray_curtis_distance_matrix.qza \
				--m-metadata-file ../mapping.tsv \
				--m-metadata-column $COMPARISON \
				--o-visualization $dir_filt_feature/bray_curtis_distance_matrix.qzv \
				--p-pairwise
			if [[ ! -f $dir_filt_feature/bray_curtis_distance_matrix.qzv ]]
			then
				touch $dir_filt_feature/bray_curtis_distance_matrix.qzv
			fi						  
			echo "DONE: diversity beta-group-significance bray_curtis_distance_matrix $dir_filt_feature $rarefaction_depth"
		fi
		dir_filt_feature=$save_dir_filt_feature_2
		if [[ ! -f $dir_filt_feature/comp-table.qza ]]
		then
			echo "RUNNING: taxa add-pseudocount"
			qiime composition add-pseudocount \
			--i-table $dir_filt_feature/table.qza \
			--o-composition-table $dir_filt_feature/comp-table.qza
			echo "DONE:  taxa add-pseudocount"
			verification_result "$dir_filt_feature/comp-table.qza"
		fi 

		if [[ ! -f $dir_filt_feature/done_diversity_vegan_TSS ]]
		then
			echo "Rscript $RSCRIPTS_DIR/diversity_vegan_tss.R $dir_filt_feature $COMPARISON $rarefaction_depth $code"
			echo "RUNNING: Rscript diversity_vegan_tss"
			conda deactivate
			Rscript $RSCRIPTS_DIR/diversity_vegan_tss.R $dir_filt_feature $COMPARISON $rarefaction_depth $code
			echo "DONE: Rscript diversity_vegan_tss"
		fi
		
		dir_filt_feature=$dir_filt_feature/$TAX_PREFIX
		if [[ ! -d  $dir_filt_feature ]]
		then
			mkdir $dir_filt_feature
		fi	
		
		if [[ ! -f $dir_filt_feature/taxa-bar-plots.qzv ]]
		then
			activate_qiime2
			echo "RUNNING:  taxa barplot $dir_filt_feature"
			qiime taxa barplot \
			--i-table $dir_filt_feature/../table.qza \
			--i-taxonomy $dir_table/$TAX_PREFIX-taxonomy.qza \
			--m-metadata-file ../mapping.tsv \
			--o-visualization $dir_filt_feature/taxa-bar-plots.qzv
			echo "DONE:  taxa barplot $dir_filt_feature"
			verification_result "$dir_filt_feature/taxa-bar-plots.qzv"
		fi
		if [[ ! -f $dir_filt_feature/taxa-bar-plots.xlsx ]]
		then
			conda deactivate
			echo "Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_filt_feature/../table.qza $dir_table/$TAX_PREFIX-taxonomy.qza $dir_filt_feature/taxa-bar-plots.xlsx"
			echo "RUNNING: Rscript extract_taxo_barplots"
			Rscript $RSCRIPTS_DIR/extract_taxo_barplots.R $dir_filt_feature/../table.qza $dir_table/$TAX_PREFIX-taxonomy.qza $dir_filt_feature/taxa-bar-plots.xlsx
			echo "DONE: Rscript extract_taxo_barplots"
			activate_qiime2
		fi	

		if [[ ! -f $dir_filt_feature/$ANCOMBC_NB_PERC/done_DifferentialAbundance ]] || [[ $BOXPLOT == "yes" && ! -f $dir_filt_feature/$ANCOMBC_NB_PERC/done_ANCOM-BC_boxplot_2 ]]
		then
			if [[ ! -d  $dir_filt_feature/$ANCOMBC_NB_PERC ]]
			then
				mkdir $dir_filt_feature/$ANCOMBC_NB_PERC
			fi	
			echo "Rscript $RSCRIPTS_DIR/run_ancombc.R $dir_filt_feature $dir_table/$TAX_PREFIX-taxonomy.qza $COMPARISON $BOXPLOT $ANCOMBC_NB_PERC"
			echo "RUNNING: Rscript ANCOM-BC"
			conda deactivate
			Rscript $RSCRIPTS_DIR/run_ancombc.R $dir_filt_feature $dir_table/$TAX_PREFIX-taxonomy.qza $COMPARISON $BOXPLOT $ANCOMBC_NB_PERC
			echo "DONE: Rscript ANCOM-BC"
		fi
	
		echo 'perl '"$PLSCRIPTS_DIR"'/merged_da_results.pl '"$dir_filt_feature"' '"$ANCOMBC_NB_PERC"
		perl $PLSCRIPTS_DIR/merged_da_results.pl $dir_filt_feature $ANCOMBC_NB_PERC;

		if [[ -f $dir_filt_feature/$ANCOMBC_NB_PERC/lefse_pairwise.tsv ]]
		then
			echo 'perl '"$PLSCRIPTS_DIR"'/merged_da_results.pl '"$dir_filt_feature"' '"$ANCOMBC_NB_PERC" 'pairwise'
			perl $PLSCRIPTS_DIR/merged_da_results.pl $dir_filt_feature $ANCOMBC_NB_PERC 'pairwise';
		fi

		echo "RUNNING: perl summary_results"
		echo 'perl '"$PLSCRIPTS_DIR"'/summary_results.pl '"$dir_filt_feature"' '"$COMPARISON"' '"$FILTER_NAME"
		perl $PLSCRIPTS_DIR/summary_results.pl $dir_filt_feature $COMPARISON $FILTER_NAME;
		echo "DONE: perl summary_results"					
	else
		echo "You must define a sampling depth at this step for $FILTER_NAME $COMPARISON $ANCOM_vs"
	fi
done
echo "-----------DONE-----------"

}