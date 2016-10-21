#!/usr/bin/perl

# Name: Position_shift_SAM.pl
# Author: Tovah Markowitz
# Date: 1/30/13
# Updated: 10/19/16 to work with telomere or rDNA pipelines
# Description: When mapping reads to a short fragment of the genome, locations must be updated to match the location in the genome
# Input: 1) SAM file with no values in chromosome location, 2) chromosome number, 3) start position of genome fragment
# Output: SAM file updated with complete information

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $help=0;

my ($inFile, $chr, $position, $headerFile, $outFile);

# set commandline options.
GetOptions ("help" => \$help,
	    "input=s" => \$inFile,
	    "chromosome=s" => \$chr,
	    "position=i" => \$position,
	    "output=s" => \$outFile) or pod2usage(2);
pod2usage(1) if $help;

################################################################
## HEADER arrays

# SAM header for SacCer3
my @SacCer3header;
$SacCer3header[0] = "\@HD\tVN:1.0\tSO:unsorted";
$SacCer3header[1] = "\@SQ\tSN:chrI\tLN:230218";
$SacCer3header[2] = "\@SQ\tSN:chrII\tLN:813184";
$SacCer3header[3] = "\@SQ\tSN:chrIII\tLN:316620";
$SacCer3header[4] = "\@SQ\tSN:chrIV\tLN:1531933";
$SacCer3header[5] = "\@SQ\tSN:chrV\tLN:576874";
$SacCer3header[6] = "\@SQ\tSN:chrVI\tLN:270161";
$SacCer3header[7] = "\@SQ\tSN:chrVII\tLN:1090940";
$SacCer3header[8] = "\@SQ\tSN:chrVIII\tLN:562643";
$SacCer3header[9] = "\@SQ\tSN:chrIX\tLN:439888";
$SacCer3header[10] = "\@SQ\tSN:chrX\tLN:745751";
$SacCer3header[11] = "\@SQ\tSN:chrXI\tLN:666816";
$SacCer3header[12] = "\@SQ\tSN:chrXII\tLN:1078177";
$SacCer3header[13] = "\@SQ\tSN:chrXIII\tLN:924431";
$SacCer3header[14] = "\@SQ\tSN:chrXIV\tLN:784333";
$SacCer3header[15] = "\@SQ\tSN:chrXV\tLN:1091291";
$SacCer3header[16] = "\@SQ\tSN:chrXVI\tLN:948066";
my $SacCer3head = join( "\n", @{SacCer3header} );

# SAM header for SK1 from the Keeney lab (2011)
my @SK1header;
$SK1header[0] = "\@HD\tVN:1.0\tSO:unsorted";
$SK1header[1] = "\@SQ\tSN:chr01\tLN:203893";
$SK1header[2] = "\@SQ\tSN:chr02\tLN:794508";
$SK1header[3] = "\@SQ\tSN:chr03\tLN:342718";
$SK1header[4] = "\@SQ\tSN:chr04\tLN:1490682";
$SK1header[5] = "\@SQ\tSN:chr05\tLN:602514";
$SK1header[6] = "\@SQ\tSN:chr06\tLN:284456";
$SK1header[7] = "\@SQ\tSN:chr07\tLN:1067526";
$SK1header[8] = "\@SQ\tSN:chr08\tLN:544538";
$SK1header[9] = "\@SQ\tSN:chr09\tLN:435585";
$SK1header[10] = "\@SQ\tSN:chr10\tLN:719294";
$SK1header[11] = "\@SQ\tSN:chr11\tLN:687260";
$SK1header[12] = "\@SQ\tSN:chr12\tLN:1008248";
$SK1header[13] = "\@SQ\tSN:chr13\tLN:908607";
$SK1header[14] = "\@SQ\tSN:chr14\tLN:812465";
$SK1header[15] = "\@SQ\tSN:chr15\tLN:1054033";
$SK1header[16] = "\@SQ\tSN:chr16\tLN:921188";
my $SK1head = join( "\n", @{SK1header} );

################################################################
## FUNCTION

# Ensures chromosome name is recognized before starting
# First determine proper header to use
my $newHead;
my @newHeader;
if ( $SK1head =~ m(SN:$chr\t) ) {
  $newHead = $SK1head;
  @newHeader = @SK1header;
} elsif ( $SacCer3head =~ m(SN:$chr\t) ) {
  $newHead = $SacCer3head;
  @newHeader = @SacCer3header;
} else {
  die("Chromosome name is not listed in the headers included with this function.\n")
}

################
## Part 1: Read in SAM file to be converted and convert
# ensure inFile exists, else die
open (IN, "<$inFile") or die ("Can not find your input file $inFile: No such file or directory.\n");

# set up necessary arrays
my @SAMData;
my @header;
my $headerCounter = -1;
my $seqCounter = -1;

# for as long as the file has more data:
while (<IN>) {
  chomp;
  # for each header line in the SAM file, save in the header array
  if ( $_ =~ /^@/ ) {
    $headerCounter++;
    $header[$headerCounter] = $_;
    # for each read, split the string into multiple locations within the array
  } else {
    $seqCounter++;
    $SAMData[$seqCounter] = [ split("\t", $_) ];
    # the chromosome number is in column 3
    $SAMData[$seqCounter][2] = $chr;
    # the start position is in column 4
    $SAMData[$seqCounter][3] = $SAMData[$seqCounter][3] + $position - 1;
  }
}
close IN;

##################
## Part 2: write output file

open OUT, ">$outFile", or die ("Can not open $outFile for writing.\n");

# if header contains @HD header line, only want to get the @SD information from this function
if ($header[0] =~ m(^\@HD) ) {
  print OUT ($header[0] . "\n");
  print OUT join( "\n", @{newHeader[1..16]} );
  print OUT "\n";
} else {
  print OUT "$newHead\n";
}
if ($header[-1] =~ m(^\@PG) ) {
  print OUT ($header[-1] . "\n");
} else {
  print ("CAREFUL: missing a line of the SAM header. Continuing anyway.\n");
}

#print the SAM output to file
for (my $x=0; $x<scalar(@SAMData);$x++) {
  my $line = join("\t",@{$SAMData[$x]});
  print OUT "$line\n";
}

close OUT;


################################################################

__END__

=head1 SYNOPSIS

Position_shift_SAM.pl [Options]

Purpose: Use this function if fastq sequences are mapped to a genome fragment
(as opposed to an extire chromosome) and the exact position of that fragment 
on the reference genome is known. 

Output: a SAM file with chromosome name included and positions shifted to reflect
the relative position of the read to the fragment and the fragment to the chromosome.

Example: rDNA in SacCer3

"perl Position_shift_SAM.pl -i IN.sam -c chrXII -p 451418 -o OUT.sam"

Options:

     -h, --help        brief help message
     -i, --input       input SAM file with no chromosome information
     -c, --chromosome  chromosome where SAM information should be mapped
     -p, --position    start position on the chromosome, shift amount
     -o, --output      name of output SAM file

=cut
