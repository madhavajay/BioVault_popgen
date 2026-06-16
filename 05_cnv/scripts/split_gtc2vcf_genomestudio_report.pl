#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my ($final_report, $out_dir, $suffix);
$suffix = '.txt';

GetOptions(
    'final-report=s' => \$final_report,
    'out-dir=s'      => \$out_dir,
    'suffix=s'       => \$suffix,
) or die "invalid arguments\n";

die "--final-report is required\n" unless defined $final_report;
die "--out-dir is required\n" unless defined $out_dir;

make_path($out_dir);

open my $in, '<', $final_report or die "cannot read $final_report: $!\n";

sub norm {
    my ($value) = @_;
    $value = lc($value // '');
    $value =~ s/^\s+|\s+$//g;
    $value =~ s{[ /]+}{_}g;
    return $value;
}

sub safe_name {
    my ($value) = @_;
    $value =~ s/^\s+|\s+$//g;
    $value =~ s{[\\/]}{_}g;
    $value =~ s/[^A-Za-z0-9_.-]+/_/g;
    return $value;
}

my @fields;
my %idx;
my $delim = "\t";

while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^\[/;
    $delim = ($line =~ /\t/) ? "\t" : ",";
    @fields = split(/\Q$delim\E/, $line, -1);
    for my $i (0 .. $#fields) {
        $idx{norm($fields[$i])} = $i;
    }
    last;
}

die "no header found in $final_report\n" unless @fields;

my $name_i = $idx{snp_name} // $idx{name};
die "header must contain Name or SNP Name\n" unless defined $name_i;

my %sample;
for my $i (0 .. $#fields) {
    my $field = $fields[$i];
    if ($field =~ /^(.+)\.B Allele Freq$/) {
        $sample{$1}{baf} = $i;
    } elsif ($field =~ /^(.+)\.Log R Ratio$/) {
        $sample{$1}{lrr} = $i;
    }
}

my @samples = sort grep {
    defined $sample{$_}{baf} && defined $sample{$_}{lrr}
} keys %sample;

die "no sample B Allele Freq / Log R Ratio column pairs found\n" unless @samples;

my %fh;
for my $sample (@samples) {
    my $safe = safe_name($sample);
    my $path = "$out_dir/$safe$suffix";
    open my $out, '>', $path or die "cannot write $path: $!\n";
    print {$out} "Name\t$sample.Log R Ratio\t$sample.B Allele Freq\n";
    $fh{$sample} = $out;
}

my $rows = 0;
my $nan_values = 0;
while (my $line = <$in>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @values = split(/\Q$delim\E/, $line, -1);
    my $name = $values[$name_i] // '';
    next unless $name ne '';
    for my $sample (@samples) {
        my $lrr = $values[$sample{$sample}{lrr}] // '';
        my $baf = $values[$sample{$sample}{baf}] // '';
        if ($lrr =~ /^-?nan$/i || $lrr eq '') {
            $lrr = 0;
            $nan_values++;
        }
        if ($baf =~ /^-?nan$/i || $baf eq '') {
            $baf = 0;
            $nan_values++;
        }
        print {$fh{$sample}} join("\t", $name, $lrr, $baf), "\n";
    }
    $rows++;
}

for my $sample (@samples) {
    close $fh{$sample};
}
close $in;

print STDERR "NOTICE: Finished processing $rows markers and generated " . scalar(@samples) . " output signal intensity files\n";
print STDERR "NOTICE: Replaced $nan_values missing/non-numeric LRR/BAF values with 0\n";
