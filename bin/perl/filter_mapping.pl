# filter_mapping.pl
# Reads filter criteria encoded as "key=value@key=value2,value3@..."
# and writes matching sample IDs from mapping.tsv to a text file.
#
# Author: Jeremy Tournayre
#
# Usage (called by qiime2.sh via eval):
#   perl filter_mapping.pl <filter_name> <filter_string>
#
# filter_string format: "col=val1,val2@col2=val3@..."
# Backslashes introduced by the orchestrator are stripped.

use strict;
use warnings;

# Parse the filter string into a hash: column => set of allowed values.
my %filter;
my @split = split('@', $ARGV[1]);
foreach (@split) {
    my $a       = $_;
    my @split2  = split("=", $a);
    my @split3  = split(",", $split2[1]);
    foreach (@split3) {
        $_ =~ s/\\//g;    # Strip backslash-escaped spaces.
        $filter{ $split2[0] }{$_} = 1;
    }
}

open(O, "> ../Filter/Samples/" . $ARGV[0] . ".txt") or die "Cannot write: $!";
print O "id\n";

open(F, "../mapping.tsv") or die "Cannot open mapping.tsv: $!";
my $header = <F>;
chomp($header);
my @title = split("\t", $header);
my $skip  = <F>;    # Skip the #q2:types row.

while (<F>) {
    chomp($_);
    my @split = split("\t", $_);
    my $i     = 0;
    my $nb    = 0;
    foreach (@split) {
        if (defined($filter{ $title[$i] })) {
            if (defined($filter{ $title[$i] }{$_})) {
                $nb++;
            }
        }
        $i++;
    }
    # A sample passes if it matches all filter criteria.
    if ($nb == keys %filter) {
        print O $split[0] . "\n";
    }
}
close(F);
close(O);