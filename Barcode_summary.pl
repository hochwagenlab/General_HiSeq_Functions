#!/usr/bin/perl

# Name: Barcode_summary.pl
# Author: Tovah Markowitz
# Date: 5/16/18

# Description: To determine the variability of barcodes within an individual
# fastq file

# Input: a fastq file
# Output: a txt file containing all barcode sequences within the fastq file
# and their frequencies

use strict;
use Data::Dumper;
use Getopt::Long;                                                               use Pod::Usage;

my $help=0;
my($fastqFile,$outputFile);

# set commandline options.
GetOptions ("fastq=s" => \$fastqFile,
	    "output=s" => \$outputFile,
	    "help" => \$help) or pod2usage(2);
pod2usage(1) if $help;

###################################

# ensure fasta file exists else die
open (IN, "<$fastqFile") or die 
    ("Cannot find your input file $fastqFile: No such file or directory\n");

# set up necessary hash
my %barcodes;
my $seqCounter=-1;

# for as long as the file has more data:
while (<IN>){
    chomp;
# place every barcode into the hash
# using the @ symbol to identify the start of each sequence
# keys are barcode sequences, values are number of times it is found
    if ($_ =~ /@/) {
	$seqCounter++;
	my $names=$_;
	my @header= split(":",$names);
	$barcodes{$header[-1]}++;
    }
}
close IN;

open (OUT, ">$outputFile") or die ("Cannot open $outputFile for writing.\n");
print OUT "Barcode\tCount\n";

foreach my $barcode (sort keys %barcodes) {
    print OUT
	"$barcode\t$barcodes{$barcode}\n";
}

close OUT;

###########################################

__END__

=head1 SYNOPSIS

Barcode_summary.pl [Options]

Purpose: To extract all the barcode sequences from a FASTQ file and
summarize to create a table of unique barcode sequences and the number of times
each has been found.

Options:

     -help/-h      brief help message
     -fastq/-f     input fastq file
     -output/-o    ouptut file name

=cut
