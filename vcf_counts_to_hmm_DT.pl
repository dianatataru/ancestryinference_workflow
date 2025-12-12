#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl "${PATH_SCRIPTS}/vcf_counts_to_hmmv3.pl" "$COUNTS_BED" "$AIM_COUNTS" 0.00000002 > "${COUNTS}.hmmsites1"

my ($counts_file, $aims_file, $rec) = @ARGV;
die "Usage: $0 <hybrid_counts_file> <AIMs_file> <recomb_rate>\n" unless $counts_file && $aims_file && defined $rec;

# Load AIMs file
my %aims;
open my $AIM, '<', $aims_file or die "Cannot open AIMs file $aims_file: $!\n";
while(<$AIM>) {
    chomp;
    my @cols = split(/\t|\s+/); # support tabs or spaces
    my $chrom = $cols[0];
    my $pos = $cols[1];
    my @panel_counts = @cols[2..$#cols];  # all remaining columns
    $aims{$chrom}{$pos} = \@panel_counts;
}
close $AIM;

# Process hybrid counts file 
my %prev_pos;  # keep track of previous position per chromosome
open my $COUNT, '<', $counts_file or die "Cannot open counts file $counts_file: $!\n";

while(<$COUNT>) {
    chomp;
    my @cols = split(/\t/);
    my $chrom = $cols[0];
    my $pos = $cols[1];

    # Skip positions not in AIMs
    next unless exists $aims{$chrom}{$pos};

    # Calculate distance
    my $dist = 0;
    if (exists $prev_pos{$chrom}) {
        $dist = ($pos - $prev_pos{$chrom}) * $rec;
    }
    $prev_pos{$chrom} = $pos;

    # Get panel counts
    my @panel_counts = @{ $aims{$chrom}{$pos} };

    # Get hybrid counts (index 4 onwards)
    my @sample_counts = @cols[4..$#cols];

    # Print final line: Chrom, Pos, panel counts, distance, sample counts
    print join("\t", $chrom, $pos, @panel_counts, $dist, @sample_counts), "\n";
}
close $COUNT;


