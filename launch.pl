#!/usr/bin/perl
#$ -S /usr/bin/perl
#$ -N JT_qiime2
#$ -V
#$ -o JT_qiime2.out
#$ -cwd
#$ -pe thread 4
#$ -q long.q
#
# launch.pl - orchestrator: reads input.tsv and runs qiime2.sh for each target.
#
# Author: Jeremy Tournayre
#
# The `#$` lines above are SGE directives; they are ignored when the script is
# run directly with `perl launch.pl` and are only used by `qsub`.
#
# Run this script from the pipeline root (where config.sh and qiime2.sh
# live) — all paths below are relative to the current working directory,
# not to this script's own location.
#
# Targets are created under Data/<qiime2_version>/<target>, where
# <qiime2_version> is CONDA_ENV_QIIME2 as set in config.sh at the moment
# this script runs.
#
# Typical usage:
#   perl launch.pl input.tsv

#### -pe thread 20
#### -q infinit.q
#### -l h_vmem=40G
#### -l mem_free=160G
use strict;
use warnings;

# Read CONDA_ENV_QIIME2 from config.sh (once, at submission time) to name
# the QIIME2-version directory targets are created under.
my $conda_env_qiime2 = `bash -c "source config.sh && echo \\\$CONDA_ENV_QIIME2"`;
chomp($conda_env_qiime2);
if (!length($conda_env_qiime2)) {
    die "Could not read CONDA_ENV_QIIME2 from config.sh. Make sure you run this script from the pipeline root (where config.sh lives).\n";
}

my $version_dir = "Data/$conda_env_qiime2";
print "QIIME2 environment (from config.sh): $conda_env_qiime2\n";
print "Targets will be created under: $version_dir/\n";

my %targets;

# Input file defaults to input.tsv; override with first argument.
my $file = "input.tsv";
if (defined($ARGV[0])) {
    $file = $ARGV[0];
}
`dos2unix $file`;

# Parse the tab-separated input file.
# Named columns populate %targets{target}{column}.
# Extra columns beyond the header are accumulated as filter key=value pairs.
open(F, $file) or die "Cannot open $file: $!";
my $row = <F>;
chomp($row);
my @title = split("\t", $row);

while (<F>) {
    chomp($_);
    my @split = split("\t", $_);
    my $target = shift(@split);
    my $i      = 1;
    my $j      = 0;
    my $title_j;
    foreach (@split) {
        $_ =~ s/ /\\ /g;    # Escape spaces in values (consumed by filter_mapping.pl).
        if (!defined($title[$i])) {
            # Columns beyond the header are filter key=value pairs.
            if ($j == 0) {
                $title_j = $_;
                $j = 1;
            } else {
                $targets{$target}{"filter"} .= $title_j . '=' . $_ . "@";
                $j = 0;
            }
        } else {
            $targets{$target}{ $title[$i] } = $_;
        }
        $i++;
    }
}
close(F);

# Required parameter columns.
my @arguments = (
    "primer_F", "primer_R", "primer_Rreverse", "primer_Freverse",
    "remove_reverse", "trim_F", "trim_R", "verif_primer_file",
    "p-min-overlap", "taxonomy_database", "taxonomy_database_primer",
    "comparison", "filter_name", "filter", "sampling_depth",
    "single_or_paired", "filtFeatures_file_nb", "filtFeatures_val",
    "boxplot", "trunc_F", "trunc_R", "rarefaction", "only_all_group",
    "p-min-fold-parent-over-abundance", "vsearch", "ANCOMBC_Nb_Perc"
);

# Maps input.tsv column names to the shell variable names used in qiime2.sh.
my %conf_map = (
    'primer_F'                         => 'PRIMER_F',
    'primer_R'                         => 'PRIMER_R',
    'primer_Rreverse'                  => 'PRIMER_RREVERSE',
    'primer_Freverse'                  => 'PRIMER_FREVERSE',
    'remove_reverse'                   => 'REMOVE_REVERSE',
    'trim_F'                           => 'TRIM_F',
    'trim_R'                           => 'TRIM_R',
    'verif_primer_file'                => 'VERIF_PRIMER_FILE',
    'p-min-overlap'                    => 'P_MIN_OVERLAP',
    'taxonomy_database'                => 'TAXONOMY_DATABASE',
    'taxonomy_database_primer'         => 'TAXONOMY_DATABASE_PRIMER',
    'comparison'                       => 'COMPARISON',
    'filter_name'                      => 'FILTER_NAME',
    'sampling_depth'                   => 'SAMPLING_DEPTH',
    'single_or_paired'                 => 'SINGLE_OR_PAIRED',
    'filtFeatures_file_nb'             => 'FILTFEATURES_MODE',
    'filtFeatures_val'                 => 'FILTFEATURES_VAL',
    'boxplot'                          => 'BOXPLOT',
    'trunc_F'                          => 'TRUNC_F',
    'trunc_R'                          => 'TRUNC_R',
    'rarefaction'                      => 'RAREFACTION',
    'only_all_group'                   => 'ONLY_ALL_GROUP',
    'p-min-fold-parent-over-abundance' => 'P_MIN_FOLD_PARENT',
    'vsearch'                          => 'VSEARCH',
    'ANCOMBC_Nb_Perc'                  => 'ANCOMBC_NB_PERC',
);

foreach (sort keys %targets) {
    my $target     = $_;
    my $target_dir = "$version_dir/$target";

    # Create output directories if this is the first run for this target.
    if (!-d $target_dir . '/QIIME2') {
        `mkdir -p $target_dir/tmp`;
        `mkdir -p $target_dir/QIIME2`;
    }
    if (!-d $target_dir . '/Filter')          { `mkdir -p $target_dir/Filter`; }
    if (!-d $target_dir . '/Filter/Samples')  { `mkdir -p $target_dir/Filter/Samples`; }
    if (!-d $target_dir . '/Filter/Features') { `mkdir -p $target_dir/Filter/Features`; }
    my $dir = $target_dir . '/QIIME2/Filter_Stat';
    if (!-d $dir) { `mkdir -p $dir`; }

    # Dedicated subdirectories for per-job run.conf and log files: several
    # input.tsv rows can share the same target and be dispatched as
    # separate concurrent qsub jobs, so each job gets its own files.
    my $logs_dir     = "$target_dir/logs";
    my $run_conf_dir = "$target_dir/run_conf";
    `mkdir -p $logs_dir`     unless -d $logs_dir;
    `mkdir -p $run_conf_dir` unless -d $run_conf_dir;

    # Validate that all required parameters are present.
    my $bool = 0;
    foreach (@arguments) {
        print $target . '--' . $_ . '--' . (defined($targets{$target}{$_}) ? $targets{$target}{$_} : '') . "\n";
        if (!defined($targets{$target}{$_})) {
            print 'Miss:' . $target . '--' . $_ . "\n";
            $bool = 1;
        }
    }

    if ($bool == 1) {
        print $target . "- ERROR: missing at least one argument\n";
        next;
    }

    # Remove the trailing '@' from the accumulated filter string.
    chop($targets{$target}{"filter"});

    # Job identifier: SGE job number (sortable by submission order/date)
    # as prefix, then the filter_name (human-readable), then a PID+random
    # suffix to guarantee uniqueness even if two rows share the same
    # filter_name. $ENV{JOB_ID} is set by SGE for qsub-submitted jobs;
    # falls back to the process PID when run directly (no qsub).
    my $sge_job_id = defined($ENV{'JOB_ID'}) ? $ENV{'JOB_ID'} : $$;
    my $job_id = $sge_job_id . '_' . $targets{$target}{'filter_name'} . '_'
               . join('', map { chr(65 + int(rand(26))) } 1..6);
    $job_id =~ s/[^-\w.]/_/g;    # sanitize for use in a filename

    my $out       = "$logs_dir/$job_id.log";
    my $run_conf  = "$run_conf_dir/$job_id.conf";

    open(my $cfh, '>', $run_conf) or die "Cannot write $run_conf: $!";
    print $cfh "# Auto-generated by launch.pl on " . localtime() . " - do not edit by hand.\n";
    foreach my $col (sort keys %conf_map) {
        my $val = defined($targets{$target}{$col}) ? $targets{$target}{$col} : '';
        $val =~ s/'/'\\''/g;    # POSIX-safe single-quote escaping.
        print $cfh $conf_map{$col} . "='" . $val . "'\n";
    }
    # FILTER is the accumulated "key=value@..." string.
    my $filt = defined($targets{$target}{"filter"}) ? $targets{$target}{"filter"} : '';
    $filt =~ s/'/'\\''/g;
    print $cfh "FILTER='" . $filt . "'\n";
    close($cfh);

    print "go\n";
    system("bash qiime2.sh '$target_dir' '$run_conf' 1>>$out 2>&1");

    # Clean up the per-job run.conf once qiime2.sh has finished with it —
    # it has already been sourced by then, and PROVENANCE.txt / the log
    # file keep a record of what ran.
    unlink($run_conf);
}