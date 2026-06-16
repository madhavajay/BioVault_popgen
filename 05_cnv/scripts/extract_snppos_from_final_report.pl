#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my ($final_report, $out);
GetOptions(
    'final-report=s' => \$final_report,
    'out=s'          => \$out,
) or die "invalid arguments\n";

die "--final-report is required\n" unless defined $final_report;
die "--out is required\n" unless defined $out;

open my $in, '<', $final_report or die "cannot read $final_report: $!\n";
open my $out_fh, '>', $out or die "cannot write $out: $!\n";

sub norm {
    my ($value) = @_;
    $value = lc($value // '');
    $value =~ s/^\s+|\s+$//g;
    $value =~ s{[ /]+}{_}g;
    return $value;
}

my $in_data = 0;
my @fields;
my %idx;
my %seen;
print {$out_fh} "Name\tChr\tPosition\n";

while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    if ($line =~ /^\[Data\]\s*$/i) {
        $in_data = 1;
        next;
    }
    next if !$in_data && $line =~ /^\[/;

    my $delim = ($line =~ /\t/) ? "\t" : ",";
    if (!@fields) {
        my @candidate = split(/\Q$delim\E/, $line, -1);
        my %candidate_idx;
        for my $i (0 .. $#candidate) {
            $candidate_idx{norm($candidate[$i])} = $i;
        }
        my $has_name = exists $candidate_idx{snp_name} || exists $candidate_idx{name};
        my $has_chr = exists $candidate_idx{chr} || exists $candidate_idx{chromosome};
        my $has_pos = exists $candidate_idx{position} || exists $candidate_idx{pos} || exists $candidate_idx{mapinfo} || exists $candidate_idx{map_info};
        next unless $has_name && $has_chr && $has_pos;
        @fields = @candidate;
        %idx = %candidate_idx;
        next;
    }

    my $name_i = $idx{snp_name} // $idx{name};
    my $chr_i = $idx{chr} // $idx{chromosome};
    my $pos_i = $idx{position} // $idx{pos} // $idx{mapinfo} // $idx{map_info};
    last unless defined $name_i && defined $chr_i && defined $pos_i;

    my @values = split(/\Q$delim\E/, $line, -1);
    my $name = $values[$name_i] // '';
    my $chr = $values[$chr_i] // '';
    my $pos = $values[$pos_i] // '';
    for my $value ($name, $chr, $pos) {
        $value =~ s/^\s+|\s+$//g;
    }
    next unless $name && $chr && $pos;
    next if $seen{$name}++;
    print {$out_fh} join("\t", $name, $chr, $pos), "\n";
}

close $in;
close $out_fh;
