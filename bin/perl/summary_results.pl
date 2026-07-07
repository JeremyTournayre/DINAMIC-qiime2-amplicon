# summary_results.pl
# Aggregates alpha and beta diversity statistics from QIIME2 core-metrics
# .qzv files and from vegan/TSS output TSVs into two summary tables:
#   results_diversity.tsv          - pairwise group comparisons
#   results_diversity_group.tsv    - per-group global tests
# plus their vegan/TSS equivalents.
#
# Author: Jeremy Tournayre
#
# Usage (called by qiime2.sh):
#   perl summary_results.pl <taxon_dir> <comparison_col> <filter_name>
#
# NOTE: ConcatVar MUST remain the last column of mapping.tsv so that
# $title[$#title] eq "ConcatVar" correctly identifies the label column.

use strict;
use warnings;

my $dir = $ARGV[0] . '/..';

# Load sample IDs for this filter.
my %list_samples;
open(F, $dir . '/sample.txt') or die "Cannot open sample.txt: $!";
while (<F>) {
    chomp($_);
    $list_samples{$_} = 1;
}
close(F);

# Build group => ConcatVar label map from mapping.tsv.
my %group2samples;
open(F, "../mapping.tsv") or die "Cannot open mapping.tsv: $!";
my $ligne = <F>;
chomp($ligne);
my @title = split("\t", $ligne);

my $col_idx = 0;
foreach (@title) {
    last if $ARGV[1] eq $_;
    $col_idx++;
}

my $skip = <F>;    # Skip the #q2:types row.
while (<F>) {
    chomp($_);
    my @split = split("\t", $_);
    next if ($ARGV[2] ne "ALL" && !defined($list_samples{ $split[0] }));

    # Use ConcatVar (last column) as the label if present.
    my $label_col = ($title[$#title] eq "ConcatVar") ? $#title : 0;
    if (!defined($group2samples{ $split[$col_idx] })) {
        $group2samples{ $split[$col_idx] } = $split[$label_col];
    } else {
        $group2samples{ $split[$col_idx] } .= ',' . $split[$label_col];
    }
}

# Helper: unzip a .qzv artifact and return its UUID code.
my ($file, $code);
sub extract_qzv {
    $file =~ s/.qzv//;
    `cp $dir/$file.qzv $dir/$file.zip`;
    `unzip -o $dir/$file.zip -d $dir/$file`;
    $code = `ls $dir/$file`;
    chomp($code);
}

# -----------------------------------------------------------------------
# QIIME2 core-metrics diversity results
# -----------------------------------------------------------------------
my %nb;    # Sample counts per group (populated from alpha significance files).

open(O2, ">" . $dir . '/results_diversity.tsv')       or die "Cannot write: $!";
print O2 "File\tGroup_1\tMean_Group_1\tSEM_Group_1\tnb_Group_1\twho_Group_1\t"
       . "Group_2\tMean_Group_2\tSEM_Group_2\tnb_Group_2\twho_Group_2\tPvalue\tQvalue\n";

open(O3, ">" . $dir . '/results_diversity_group.tsv') or die "Cannot write: $!";
print O3 "File\tGroup\tnb\tPvalue\n";

# Alpha diversity: group-significance .qzv files.
my @ls2 = `ls $dir/core-metrics-results_*/*group-signifi*qzv`;
foreach (@ls2) {
    chomp($_);
    $file = $_;
    my $dir_regexp = $dir;
    $dir_regexp =~ s/\+/\\\+/g;
    $file =~ s/$dir_regexp\///;
    next if -z $dir . '/' . $file;
    extract_qzv();

    if (-e $dir . '/' . $file . '/' . $code . '/data/column-' . $ARGV[1] . '.jsonp') {
        open(F, $dir . '/' . $file . '/' . $code . '/data/column-' . $ARGV[1] . '.jsonp');
        my $line = <F>;

        # Parse group names and sample values from the jsonp file.
        $line =~ /"index":\[(.*?)\]/;
        my @groups = map { s/"([^"]+) \(\s*n=\d+\s*\)"/$1/; $_ } split(/,/, $1);
        $line =~ /"data":\[(.*?)\]\]/;
        my @data_blocks = split(/\],\[/, $1);

        my %stats;
        for my $i (0 .. $#groups) {
            my @values = $data_blocks[$i] =~ /([\d\.]+)/g;
            my ($sum, $sum_sq, $count) = (0, 0, scalar(@values));
            $sum    += $_ for @values;
            $sum_sq += $_ ** 2 for @values;
            my $mean     = $count ? $sum / $count : 'NA';
            my $variance = ($count > 1) ? ($sum_sq / $count - $mean ** 2) : 'NA';
            my $std_dev  = ($variance ne 'NA') ? sqrt($variance) : 'NA';
            my $sem      = ($count > 1) ? $std_dev / sqrt($count) : 'NA';
            $stats{ $groups[$i] } = { mean => $mean, sem => $sem };
        }
        close(F);

        open(F, $dir . '/' . $file . '/' . $code . '/data/kruskal-wallis-pairwise-' . $ARGV[1] . '.csv');
        my $skip2 = <F>;
        while (<F>) {
            chomp($_);
            my @split = split(',', $_);
            $split[0] =~ m/(.*)\s\(n=(.*)\)/;
            my ($g1, $n1) = ($1, $2);
            print O2 $file . "\t" . $g1 . "\t" . $stats{$g1}{mean} . "\t" . $stats{$g1}{sem}
                          . "\t" . $n1 . "\t" . $group2samples{$g1} . "\t";
            $nb{$g1} = $n1;

            $split[1] =~ m/(.*)\s\(n=(.*)\)/;
            my ($g2, $n2) = ($1, $2);
            $nb{$g2} = $n2;
            print O2 $g2 . "\t" . $stats{$g2}{mean} . "\t" . $stats{$g2}{sem}
                         . "\t" . $n2 . "\t" . $group2samples{$g2}
                         . "\t" . $split[3] . "\t" . $split[4] . "\n";
        }
        close(F);

        open(F, $dir . '/' . $file . '/' . $code . '/data/column-' . $ARGV[1] . '.jsonp');
        while (<F>) {
            chomp($_);
            $_ =~ m/"filtered": (.*?)\}/;
            print O3 $file . "\t" . $ARGV[1] . "\t" . $1 . "\t";
            $_ =~ m/"p": (.*?)\}/;
            print O3 $1 . "\n";
        }
        close(F);
        `rm -R $dir/$file`;
        `rm $dir/$file.zip`;
    }
}

# Beta diversity: distance-matrix group-significance .qzv files.
@ls2 = `ls $dir/core-metrics-results_*/*distance_matrix*qzv`;
foreach (@ls2) {
    chomp($_);
    $file = $_;
    my $dir_regexp = $dir;
    $dir_regexp =~ s/\+/\\\+/g;
    $file =~ s/$dir_regexp\///;
    next if -z $dir . '/' . $file;
    extract_qzv();

    if (-e $dir . '/' . $file . '/' . $code . '/data/permanova-pairwise.csv') {
        open(F, $dir . '/' . $file . '/' . $code . '/data/permanova-pairwise.csv');
        my $skip2 = <F>;
        # permanova-pairwise.csv columns: Group 1, Group 2, Sample size,
        # Permutations, pseudo-F, p-value, q-value  (indices 0..6)
        while (<F>) {
            chomp($_);
            my @split = split(',', $_);
            $split[0] =~ s/\s\(.*//;
            $split[1] =~ s/\s\(.*//;
            print O2 $file . "\t" . $split[0] . "\tBeta\tBeta\t" . $nb{ $split[0] } . "\t" . $group2samples{ $split[0] }
                         . "\t" . $split[1] . "\tBeta\tBeta\t" . $nb{ $split[1] } . "\t" . $group2samples{ $split[1] }
                         . "\t" . $split[5] . "\t" . $split[6] . "\n";
        }
        close(F);

        open(F, $dir . '/' . $file . '/' . $code . '/data/index.html');
        print O3 $file . "\t" . $ARGV[1] . "\t";
        my ($bool_0, $bool_00) = (0, 0);
        while (<F>) {
            chomp($_);
            if ($bool_0 == 1) { $_ =~ m/<td>(.*?)<\/td>/; print O3 "\t" . $1 . "\n"; $bool_0 = 2; }
            if ($bool_00)     { $_ =~ m/<td>(.*?)<\/td>/; print O3 $1; $bool_00 = 0; }
            if ($_ =~ m/<th>sample size<\/th>/) { $bool_00 = 1; }
            if ($bool_0 == 0 && $_ =~ m/<th>p-value<\/th>/) { $bool_0 = 1; }
        }
        close(F);
        `rm -R $dir/$file`;
        `rm $dir/$file.zip`;
    }
}
close(O2);
close(O3);

# -----------------------------------------------------------------------
# Vegan / TSS diversity results
# -----------------------------------------------------------------------
open(O2, ">" . $dir . '/results_diversity_vegan_TSS.tsv') or die "Cannot write: $!";
print O2 "File\tGroup_1\tMean_Group_1\tSEM_Group_1\tnb_Group_1\twho_Group_1\t"
       . "Group_2\tMean_Group_2\tSEM_Group_2\tnb_Group_2\twho_Group_2\tPvalue\tQvalue\n";

open(O3, ">" . $dir . '/results_diversity_group_vegan_TSS.tsv') or die "Cannot write: $!";
print O3 "File\tGroup\tnb\tPvalue\n";

my $nsamples = scalar keys %list_samples;

# PERMANOVA results from vegan.
foreach (`ls $dir/diversity_vegan_TSS/permanova_results*`) {
    chomp($_);
    $file = $_;
    my $dr = $dir; $dr =~ s/\+/\\\+/g;
    $file =~ s/$dr\///;
    open(F, $dir . '/' . $file) or next;
    my $skip2 = <F>;
    while (<F>) {
        chomp($_);
        if ($_ =~ m/Model/) {
            my @split = split('\t', $_);
            print O3 $split[8] . '-' . $split[7] . "\t" . $ARGV[1] . "\t" . $nsamples . "\t" . $split[5] . "\n";
        }
    }
    close(F);
}

# Kruskal-Wallis results from vegan.
foreach (`ls $dir/diversity_vegan_TSS/kruskal*`) {
    chomp($_);
    $file = $_;
    my $dr = $dir; $dr =~ s/\+/\\\+/g;
    $file =~ s/$dr\///;
    open(F, $dir . '/' . $file) or next;
    my $skip2 = <F>;
    while (<F>) {
        chomp($_);
        my @split = split('\t', $_);
        print O3 $split[7] . '-' . $split[2] . "\t" . $ARGV[1] . "\t" . $nsamples . "\t" . $split[4] . "\n";
    }
    close(F);
}
close(O3);

# Pairwise PERMANOVA (beta diversity).
foreach (`ls $dir/diversity_vegan_TSS/pairwise_results*`) {
    chomp($_);
    $file = $_;
    my $dr = $dir; $dr =~ s/\+/\\\+/g;
    $file =~ s/$dr\///;
    open(F, $dir . '/' . $file) or next;
    my $skip2 = <F>;
    while (<F>) {
        chomp($_);
        my @split  = split('\t', $_);
        my @split2 = split(' vs ', $split[1]);
        print O2 $split[12] . '-' . $split[11]
                . "\t" . $split2[0] . "\tBeta\tBeta\t" . $nb{ $split2[0] } . "\t" . $group2samples{ $split2[0] }
                . "\t" . $split2[1] . "\tBeta\tBeta\t" . $nb{ $split2[1] } . "\t" . $group2samples{ $split2[1] }
                . "\t" . $split[6] . "\t" . $split[9] . "\n";
    }
    close(F);
}

# Pairwise Wilcoxon (alpha diversity).
foreach (`ls $dir/diversity_vegan_TSS/pairwise_wilcox_results*`) {
    chomp($_);
    $file = $_;
    my $dr = $dir; $dr =~ s/\+/\\\+/g;
    $file =~ s/$dr\///;

    # Load corresponding mean values.
    my %means;
    (my $file2 = $file) =~ s/pairwise_wilcox_results/mean_values/;
    open(F, $dir . '/' . $file2) or next;
    my $skip2 = <F>;
    while (<F>) {
        chomp($_);
        my @split = split('\t', $_);
        $means{ $split[3] }{ $split[1] } = $split[2];
    }
    close(F);

    # Load corresponding SEM values.
    my %sem;
    (my $file3 = $file) =~ s/pairwise_wilcox_results/sem_values/;
    open(F, $dir . '/' . $file3) or next;
    $skip2 = <F>;
    while (<F>) {
        chomp($_);
        my @split = split('\t', $_);
        $sem{ $split[3] }{ $split[1] } = $split[2];
    }
    close(F);

    open(F, $dir . '/' . $file) or next;
    my $title = <F>;
    my @spl_title = split('\t', $title);
    shift(@spl_title); shift(@spl_title); shift(@spl_title);
    pop(@spl_title);

    while (<F>) {
        chomp($_);
        my @split = split('\t', $_);
        $split[0] =~ /(.*)\.([^.]*)$/;
        my $group = $2;
        my $i2 = 3;
        foreach (@spl_title) {
            my $a_title = $_;
            $a_title =~ s/p\.values\.//g;
            $a_title =~ s/\./ /g;
            $a_title =~ s/ (\d)/-$1/g;
            # If the R rowname had no disambiguating suffix (e.g. each group
            # appeared only once, so R never appended ".1"/".2"), the regex
            # above does not match and $group stays undef. Fall back to
            # whichever known group is NOT the current column's group.
            if (!defined($group)) {
                foreach (keys %group2samples) {
                    if ($_ eq $a_title) {
                        next;
                    } else {
                        $group = $_;
                    }
                }
            }
            next if $split[$i2] eq "NA";
            print O2 $split[-1] . '-' . $split[2]
                    . "\t" . $group     . "\t" . $means{ $split[2] }{$group}   . "\t" . $sem{ $split[2] }{$group}
                    . "\t" . $nb{$group} . "\t" . $group2samples{$group}
                    . "\t" . $a_title   . "\t" . $means{ $split[2] }{$a_title} . "\t" . $sem{ $split[2] }{$group}
                    . "\t" . $nb{$a_title} . "\t" . $group2samples{$a_title}
                    . "\tNA\t" . $split[$i2] . "\n";
            $i2++;
        }
    }
    close(F);
}
close(O2);