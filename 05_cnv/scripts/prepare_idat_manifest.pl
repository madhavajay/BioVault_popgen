#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Spec;

my ($idat_dir, $sample_map, $out, $missing_out);
GetOptions(
    'idat-dir=s'    => \$idat_dir,
    'sample-map=s'  => \$sample_map,
    'out=s'         => \$out,
    'missing-out=s' => \$missing_out,
) or die "invalid arguments\n";

die "--idat-dir is required\n" unless defined $idat_dir;
die "--out is required\n" unless defined $out;
die "--missing-out is required\n" unless defined $missing_out;

sub norm_header {
    my ($value) = @_;
    $value = lc($value // '');
    $value =~ s/^\s+|\s+$//g;
    $value =~ s{[ /]+}{_}g;
    return $value;
}

sub read_sample_map {
    my ($path) = @_;
    my %rows;
    return %rows unless defined $path && -f $path;

    open my $fh, '<', $path or die "cannot read sample map $path: $!\n";
    my $header = <$fh>;
    return %rows unless defined $header;
    chomp $header;
    my $delim = ($header =~ /\t/) ? "\t" : ",";
    my @fields = map { norm_header($_) } split(/\Q$delim\E/, $header, -1);
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @values = split(/\Q$delim\E/, $line, -1);
        my %row;
        @row{@fields} = @values;
        my $chip_id =
            ($row{chip_id_idat} // '') ||
            ($row{chip_id} // '') ||
            ($row{sentrix_id} // '') ||
            ($row{idat} // '');
        $chip_id =~ s/^\s+|\s+$//g;
        next unless $chip_id;
        my $sample_id =
            ($row{sample_id} // '') ||
            ($row{sample} // '') ||
            ($row{sample_name} // '') ||
            $chip_id;
        my $population =
            ($row{chechen_circassian} // '') ||
            ($row{population} // '') ||
            ($row{cohort} // '');
        my $sex = ($row{sex} // '') || ($row{gender} // '');
        for my $value ($sample_id, $population, $sex) {
            $value =~ s/^\s+|\s+$//g;
        }
        $rows{$chip_id} = {
            sample_id  => $sample_id,
            population => $population,
            sex        => $sex,
        };
    }
    close $fh;
    return %rows;
}

opendir my $dh, $idat_dir or die "cannot read IDAT dir $idat_dir: $!\n";
my %pairs;
while (my $name = readdir $dh) {
    next unless $name =~ /\.idat$/;
    if ($name =~ /^(.*)_Red\.idat$/) {
        $pairs{$1}{red} = File::Spec->catfile($idat_dir, $name);
    } elsif ($name =~ /^(.*)_Grn\.idat$/) {
        $pairs{$1}{green} = File::Spec->catfile($idat_dir, $name);
    }
}
closedir $dh;

my %sample_meta = read_sample_map($sample_map);

open my $out_fh, '>', $out or die "cannot write $out: $!\n";
open my $miss_fh, '>', $missing_out or die "cannot write $missing_out: $!\n";

print {$out_fh} join("\t", qw(chip_id_idat sample_id population sex red_idat green_idat)), "\n";
print {$miss_fh} join("\t", qw(chip_id_idat problem red_idat green_idat)), "\n";

for my $chip_id (sort keys %pairs) {
    my $red = $pairs{$chip_id}{red} // '';
    my $green = $pairs{$chip_id}{green} // '';
    if (!$red || !$green) {
        my $problem = $red ? 'missing_green' : 'missing_red';
        print {$miss_fh} join("\t", $chip_id, $problem, $red, $green), "\n";
        next;
    }
    my $meta = $sample_meta{$chip_id} // {};
    print {$out_fh} join(
        "\t",
        $chip_id,
        $meta->{sample_id} // $chip_id,
        $meta->{population} // '',
        $meta->{sex} // '',
        $red,
        $green,
    ), "\n";
}

close $out_fh;
close $miss_fh;
