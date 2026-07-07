#!/usr/bin/perl
# db/launch_db_multipleqsub.pl
# Splits db/input.tsv into one single-row file per classifier target, then
# submits one qsub job per target so each classifier build runs in parallel.
#
# Author: Jeremy Tournayre
#
# IMPORTANT: activate the QIIME2 environment you want to build under BEFORE
# running this script:
#   source ../config.sh
#   activate_qiime2
#   perl launch_db_multipleqsub.pl input.tsv
#
# Each qsub job starts in a fresh shell that does NOT inherit CONDA_ENV_QIIME2
# from your interactive shell. This script reads it here (from the currently
# active environment) and passes it explicitly to every submitted job via
# `qsub -v`, so all jobs build under the SAME QIIME2 version you activated.
#
# Usage:  perl launch_db_multipleqsub.pl [input_file]   (defaults to input.tsv)
#
# Each submitted job runs: perl launch_db.pl input_arg_qsub/<input_file>/<N>
# which builds exactly one classifier, using the same run.conf-based
# interface as a normal launch_db.pl call.

use strict;
use warnings;

# Capture the QIIME2 environment the person activated in THIS shell, before
# any qsub job is submitted. Each job will be told to use this exact value.
my $conda_env_qiime2 = $ENV{'CONDA_ENV_QIIME2'};
if (!defined($conda_env_qiime2) || !length($conda_env_qiime2)) {
    die "CONDA_ENV_QIIME2 is not set in this shell. Run 'source ../config.sh && activate_qiime2' first.\n";
}
print "Submitting all jobs with QIIME2 environment: $conda_env_qiime2\n";

my $file = "input.tsv";
if (defined($ARGV[0])) { $file = $ARGV[0]; }

open(F, $file) or die "Cannot open $file: $!";
my $header = <F>;

my $split_dir = "input_arg_qsub/" . $file;
system("mkdir -p $split_dir");

my $i = 0;
while (<F>) {
    chomp($_);
    $i++;
    open(O, ">" . $split_dir . "/" . $i) or die "Cannot write split file: $!";
    print O $header;
    print O $_ . "\n";
    close(O);

    # Submit one job per DB target; each job runs launch_db.pl on its
    # single-row file. -v CONDA_ENV_QIIME2=... explicitly injects the
    # captured value into the job's environment, so it does NOT fall back
    # to config.sh's default inside the fresh qsub shell.
    system("qsub -v CONDA_ENV_QIIME2=$conda_env_qiime2 launch_db.pl " . $split_dir . "/" . $i);
}
close(F);