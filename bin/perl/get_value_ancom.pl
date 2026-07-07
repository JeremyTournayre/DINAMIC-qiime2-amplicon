# get_value_ancom.pl
# Reads the sample list for a given filter, looks up the comparison column in
# mapping.tsv, and writes per-pair sample files used by the diversity steps.
# Prints the list of comparisons (e.g. "A_vs_B C_vs_D All") to stdout,
# which qiime2.sh captures to iterate over pairwise comparisons.
#
# Author: Jeremy Tournayre
#
# Usage (called by qiime2.sh):
#   perl get_value_ancom.pl <filter_name> <comparison_col> <out_dir> <only_all>
#
#   only_all: "no" => generate pairwise directories; anything else => skip pairs.

use strict;
use warnings;

# Load the pre-filtered sample IDs.
my %samples;
open(F, "../Filter/Samples/" . $ARGV[0] . ".txt") or die "Cannot open sample list: $!";
while (<F>) {
    chomp($_);
    $samples{$_} = 1;
}
close(F);

`dos2unix ../mapping.tsv`;
open(F, "../mapping.tsv") or die "Cannot open mapping.tsv: $!";
my $row = <F>;
chomp($row);
my @title = split("\t", $row);

# Find the column index of the comparison variable.
my $i = 0;
foreach (@title) {
    last if $_ eq $ARGV[1];
    $i++;
}

# Build a map: group_value => concatenated sample IDs.
my %values;
while (<F>) {
    chomp($_);
    my @split = split("\t", $_);
    next unless defined($samples{ $split[0] });
    $values{ $split[$i] } .= $split[0] . "\n";
}
close(F);

my $res = "";

if ($ARGV[3] eq "no") {
    # Generate one directory and sample file per pair.
    foreach my $a (sort keys %values) {
        foreach my $b (sort keys %values) {
            if ($a lt $b) {
                $res .= $a . '_vs_' . $b . " ";
                my $dir = $ARGV[2] . '/' . $a . '_vs_' . $b;
                `mkdir $dir` unless -d $dir;

                # Combined sample file (used by feature-table filter-samples).
                open(O, ">" . $dir . "/samples.txt") or die "Cannot write: $!";
                my $tmp_b = $values{$b};
                chop($tmp_b);
                print O "id\n";
                print O $values{$a} . $tmp_b;
                close(O);

                # Individual group files (used by diversity and DA scripts).
                open(O, ">" . $dir . "/samples_groupA.txt") or die "Cannot write: $!";
                my $tmp_a = $values{$a};
                chop($tmp_a);
                print O $tmp_a;
                close(O);

                open(O, ">" . $dir . "/samples_groupB.txt") or die "Cannot write: $!";
                my $tmp_b2 = $values{$b};
                chop($tmp_b2);
                print O $tmp_b2;
                close(O);
            }
        }
    }
}

# "All" is always appended last (load-bearing: qiime2.sh processes it last
# so that the conda environment state at loop end is predictable).
print $res . " All";