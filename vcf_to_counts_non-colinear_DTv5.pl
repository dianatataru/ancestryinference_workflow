 #!/usr/bin/perl
use strict;
use warnings;

# usage:
#   perl vcf_to_counts_non-colinear_DTv5.pl infile_aims outfile.counts

my ($infile, $outfile) = @ARGV;
die "Usage: $0 input.vcf.aims output.counts\n" unless $infile && $outfile;

open(my $IN,  "<", $infile)  or die "Cannot open $infile: $!";
open(my $OUT, ">", $outfile) or die "Cannot write $outfile: $!";

while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ /^#/;     # skip header lines
    next if $line =~ /^\s*$/;  # skip blank

    my @f = split(/\t/, $line);

    # Expected columns:
    # 0 compound
    # 1 chromosome
    # 2 position
    # 3 AIMS_REF
    # 4 AIMS_ALT
    # 5 compound (again)
    # 6 vcf_chrom
    # 7 vcf_pos
    # 8 .
    # 9 vcf_REF
    # 10 vcf_ALT
    # 11 QUAL
    # 12 FILTER
    # 13 INFO
    # 14 sample1 GT:PL
    # 15 sample2 GT:PL
    # ... etc

    next unless @f >= 15;

    my $compound  = $f[0];
    my $aims_ref  = $f[3];
    my $aims_alt  = $f[4];

    my $vcf_ref   = $f[9];
    my $vcf_alt   = $f[10];

    # Make sure REF matches
    next if $vcf_ref ne $aims_ref;

    # Make sure AIMs_ALT is present in VCF_ALT list
    my @alts = $vcf_alt eq '.' ? () : split(/,/, $vcf_alt);
    my %alt_index;
    for my $i (0 .. $#alts) {  # index 0 → ALT1 (allele "1")
        $alt_index{$alts[$i]} = $i + 1;
    }

    next unless exists $alt_index{$aims_alt};  # target ALT not present

    my $target_allele = $alt_index{$aims_alt};  # 1-based allele index

    # Print leading columns
    print $OUT join("\t", $compound, $aims_ref, $aims_alt);

    # Process each sample 
    for my $i (14 .. $#f) {
        my $sample = $f[$i];
        my ($gt) = split(/:/, $sample);  # extract genotype string

        my ($ref_count, $alt_count) = (0, 0);

        if ($gt ne "./." && $gt ne "." && $gt ne "./") {
            $gt =~ s/\|/\//g;
            my @alleles = split(/\//, $gt);

            # diploid expected; if not, treat as missing
            if (@alleles == 2) {
                for my $a (@alleles) {
                    if ($a eq ".") {
                        # missing allele → ignore
                    }
                    elsif ($a == 0) {
                        $ref_count++;
                    }
                    elsif ($a == $target_allele) {
                        $alt_count++;
                    }
                    else {
                        # allele corresponds to some OTHER ALT allele → ignore it
                    }
                }
            }
        }

        print $OUT "\t$ref_count\t$alt_count";
    }

    print $OUT "\n";
}

close $IN;
close $OUT;

print "Finished writing counts to $outfile\n";

