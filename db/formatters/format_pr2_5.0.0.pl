#!/usr/bin/perl
# db/formatters/format_pr2_5.0.0.pl
# Converts PR2 v5.0.0 UTAX-formatted FASTA into:
#   - formatted_pr2_version_5.0.0_SSU_UTAX.fasta  (clean FASTA, gene IDs only)
#   - taxonomy_pr2_version_5.0.0_SSU_UTAX.fasta   (tab-sep ID <TAB> taxonomy)
#
# Author: Jeremy Tournayre
#
# Input: pr2_version_5.0.0/pr2_version_5.0.0_SSU_UTAX.fasta
# Run this script once manually before format_pr2_5.0.0.sh.

use strict;
use warnings;

open(F,  "pr2_version_5.0.0/pr2_version_5.0.0_SSU_UTAX.fasta")          or die $!;
open(O,  ">pr2_version_5.0.0/formatted_pr2_version_5.0.0_SSU_UTAX.fasta") or die $!;
open(O2, ">pr2_version_5.0.0/taxonomy_pr2_version_5.0.0_SSU_UTAX.fasta")  or die $!;

while (<F>) {
    chomp($_);
    if ($_ =~ m/^>(.*?);tax=(.*)/) {
        my $gene_id = $1;
        my $taxo    = $2;
        # Convert UTAX separators to QIIME2 format.
        $taxo =~ s/:/__/g;
        $taxo =~ s/,/;/g;
        print O2 $gene_id . "\t" . $taxo . "\n";
        print O  '>' . $gene_id . "\n";
    } else {
        print O $_ . "\n";
    }
}
close(F);
close(O);
close(O2);