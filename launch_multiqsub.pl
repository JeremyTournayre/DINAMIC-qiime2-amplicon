#!/usr/bin/perl
# launch_multipleqsub.pl
# Splits the input file into one single-row file per target, then submits
# one qsub job per target so each analysis runs in parallel on the cluster.
#
# Author: Jeremy Tournayre
#
# Usage:  perl launch_multipleqsub.pl [input_file]   (defaults to input.txt)
#
# Each submitted job runs: perl launch.pl input_arg_qsub/<input_file>/<N>
# which processes exactly one target, using the same run.conf-based interface
# as a normal launch.pl call.

use strict;
use warnings;
use FindBin qw($Bin);

# Capture the QIIME2 environment the person activated in THIS shell, before
# any qsub job is submitted, rather than re-reading config.sh independently
# (which could have changed in the meantime).
my $conda_env_qiime2 = $ENV{'CONDA_ENV_QIIME2'};
if (!defined($conda_env_qiime2) || !length($conda_env_qiime2)) {
    die "CONDA_ENV_QIIME2 is not set in this shell. Run 'source config.sh && activate_qiime2' first.\n";
}
print "QIIME2 environment (from shell): $conda_env_qiime2\n";

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

    # Submit one job per target; each job runs launch.pl on its single-row
    # file. launch.pl is referenced via its absolute path ($Bin = the
    # directory containing THIS script, i.e. the pipeline root) so the job
    # finds it regardless of the current working directory. Each submitted
    # launch.pl independently reads config.sh to determine CONDA_ENV_QIIME2.
    system("qsub $Bin/launch.pl " . $split_dir . "/" . $i);
}
close(F);