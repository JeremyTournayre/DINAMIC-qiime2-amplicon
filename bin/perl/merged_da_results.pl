#!/usr/bin/perl
# merged_da_results.pl
# Merges differential abundance results from ANCOM-BC, DESeq2, IndicSpecies
# and LEfSe (optional) into a single TSV file.
#
# Author: Jeremy Tournayre
#
# Usage (called by qiime2.sh):
#   perl merged_da_results.pl <base_dir> <sub_dir>              # global results
#   perl merged_da_results.pl <base_dir> <sub_dir> pairwise     # pairwise results

use strict;
use warnings;
use Text::CSV;
use File::Find;
use Cwd 'abs_path';

my $pairwise = "";
if (defined($ARGV[2])) { $pairwise = 1; }

my $base_dir = shift @ARGV;
my $sub_dir  = shift @ARGV;
if ($pairwise ne "") { $pairwise = '_pairwise'; }

my $dir = "$base_dir/$sub_dir";

# Locate the ANCOM-BC output directory (named ANCOM_BC_ref_is_<group>).
my $ancom_bc_dir;
my $reference_group;
find(
    sub {
        if (-d $_ && $_ =~ /^ANCOM_BC_ref_is/) {
            $ancom_bc_dir = $File::Find::name;
            if ($_ =~ /ANCOM_BC_ref_is_([-\w]+)/) {
                $reference_group = $1;
            }
        }
    },
    $dir
);
die "Error: cannot find ANCOM_BC_ref_is* directory in $dir\n" unless $ancom_bc_dir;

my %file_paths = (
    "DESEQ2.tsv"       => "$dir/DESEQ2"       . $pairwise . ".tsv",
    "ANCOM-BC.tsv"     => "$ancom_bc_dir/ANCOM-BC" . $pairwise . ".tsv",
    "IndicSpecies.tsv" => "$dir/IndicSpecies"  . $pairwise . ".tsv",
    "lefse.tsv"        => "$dir/lefse"         . $pairwise . ".tsv",
);

# LEfSe is optional.
delete $file_paths{"lefse.tsv"} unless -e $file_paths{"lefse.tsv"};

my %data;
my @all_columns = ("name", "Stat on taxa");

my %list_groups;
my %mean_values;

if ($pairwise ne "") {
    # Collect group names from DESeq2 or IndicSpecies output.
    for my $src ("DESEQ2.tsv", "IndicSpecies.tsv") {
        next unless -e $file_paths{$src};
        open my $fh, "<", $file_paths{$src} or die "Cannot open $file_paths{$src}: $!";
        my $csv = Text::CSV->new({ sep_char => "\t", binary => 1, auto_diag => 1 });
        $csv->getline($fh);
        while (my $row = $csv->getline($fh)) {
            last unless defined($row->[2]);
            $list_groups{ $row->[2] } = 1;
            $list_groups{ $row->[3] } = 1;
        }
        last if keys %list_groups;
    }

    # Collect mean abundance column positions from ANCOM-BC output.
    open my $fh, "<", $file_paths{"ANCOM-BC.tsv"} or die "Cannot open ANCOM-BC.tsv: $!";
    my $csv = Text::CSV->new({ sep_char => "\t", binary => 1, auto_diag => 1 });
    my @cols = @{ $csv->getline($fh) };
    for my $col (0 .. $#cols) {
        if ($cols[$col] =~ /^Mean_([-\w]+)$/) {
            $mean_values{$1} = $col;
        }
    }
    close $fh;
}

my %already_done;

foreach my $file (sort keys %file_paths) {
    my $file_path = $file_paths{$file};
    open my $fh, "<", $file_path or die "Cannot open $file_path: $!";
    my $csv     = Text::CSV->new({ sep_char => "\t", binary => 1, auto_diag => 1 });
    my @columns = @{ $csv->getline($fh) };
    next unless defined($columns[1]);

    my ($group1_col, $group2_col);
    if ($file eq "DESEQ2.tsv" || $file eq "IndicSpecies.tsv" || $file eq "lefse.tsv") {
        ($group1_col, $group2_col) = (2, 3);
    }

    my @data_cols = grep { $_ > 1 } (0 .. $#columns);

    unless ($pairwise) {
        foreach my $i (@data_cols) {
            my $col_name = "${file}_$columns[$i]";
            push @all_columns, $col_name unless grep { $_ eq $col_name } @all_columns;
        }
    }

    while (my $row = $csv->getline($fh)) {
        my $name         = $row->[0];
        my $stat_on_taxa = $row->[1];

        if ($pairwise && $file eq "ANCOM-BC.tsv") {
            foreach my $col (@data_cols) {
                next if $columns[$col] =~ /^Mean_([-\w]+)$/;
                my @g;
                my $nb_group  = 0;
                my $unique_key = $name;
                foreach (sort keys %list_groups) {
                    if ($columns[$col] =~ /$_/) {
                        $g[$nb_group] = $_;
                        $unique_key  .= "_$_";
                        $nb_group++;
                    }
                }
                if ($nb_group == 1) {
                    $g[1]        = $reference_group;
                    $unique_key .= "_$reference_group";
                    $nb_group++;
                }
                next if $nb_group != 2;

                $data{$unique_key}{"name"}        = $name;
                $data{$unique_key}{"Stat on taxa"} = $stat_on_taxa;

                for my $tag ("group1", "group2") {
                    my $cn = "${file}_${tag}";
                    $data{$unique_key}{$cn} = $tag eq "group1" ? $g[0] : $g[1];
                    push @all_columns, $cn unless grep { $_ eq $cn } @all_columns;
                }
                for my $tag ("Mean_group1", "Mean_group2") {
                    my $cn  = "${file}_${tag}";
                    my $grp = $tag eq "Mean_group1" ? $g[0] : $g[1];
                    $data{$unique_key}{$cn} = $row->[ $mean_values{$grp} ];
                    push @all_columns, $cn unless grep { $_ eq $cn } @all_columns;
                }

                my $tmp = $columns[$col];
                foreach (@g) { $tmp =~ s/(.*)_.*$_/$1/ }
                my $cn = "${file}_$tmp";
                $already_done{$cn} = 1;
                push @all_columns, $cn unless grep { $_ eq $cn } @all_columns;
                $data{$unique_key}{$cn} = $row->[$col];
            }
        } else {
            my $key;
            if ($pairwise) {
                $key = "${name}_" . $row->[$group1_col] . '_' . $row->[$group2_col];
                unless (defined($data{$key})) {
                    $key = "${name}_" . $row->[$group2_col] . '_' . $row->[$group1_col];
                }
            } else {
                $key = $name;
            }
            $data{$key}{"name"}        = $name;
            $data{$key}{"Stat on taxa"} = $stat_on_taxa;
            foreach my $col (@data_cols) {
                my $cn = "${file}_$columns[$col]";
                push @all_columns, $cn unless grep { $_ eq $cn } @all_columns;
                $already_done{$cn}  = 1;
                $data{$key}{$cn}    = $row->[$col];
            }
        }
    }
    close $fh;
}

my $output_file = $pairwise ? "$dir/merged_results_pairwise.tsv" : "$dir/merged_results.tsv";
open my $out, ">", $output_file or die "Cannot write $output_file: $!";
print $out join("\t", @all_columns), "\n";

foreach my $key (sort keys %data) {
    my @row;
    foreach my $col (@all_columns) {
        push @row, $data{$key}{$col} // "NA";
    }
    print $out join("\t", @row), "\n";
}
close $out;

print "Merge complete: '$output_file' created.\n";