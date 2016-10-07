#!/usr/bin/perl

# New name: Fastq_Read_Extract_byReadName.pl
# Based upon: fastqConverterHash2.pl
# Author: Tovah Markowitz
# Date: 10/14/14

# Description: extract specific queries from a fastq file
# Input: name of fastq file, query file, output file name
# Output: query in fastQ format

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $help=0;
my ($fastqFile, $queryFile, $outputFile);

# set commandline options
GetOptions ("input=s" => \$fastqFile,
	    "query=s" => \$queryFile,
	    "output=s" => \$outputFile,,
	    "help" => \$help,) or pod2usage(2);
pod2usage(1) if $help;

#############################################

my $sequences=LoadFastq($fastqFile);

# open query file and read in queries, save only first column
# if only one column, everything becomes a query
open (QUERY, "<$queryFile") or die ("Can not find your query file $queryFile: No such file or directory\n");

my $seqCounter = -1;
my @queries;
my $query;

# read in queries
while (<QUERY>){
    chomp;
    $seqCounter++;
    $query= $_;
    $query =~ /^\S+/;
#print Dumper ($&);
    $queries[$seqCounter] = $&;
}
close QUERY;

# open output file, search and write concurrently
open (OUT, ">$outputFile") or die ("Can not open $outputFile for writing.\n"); 

# look through each query
# ensure that query is in the fastq file
# write output
for (my $q=0; $q < scalar(@queries); $q++) {
    my $fastaName = $queries[$q];
    if ($sequences -> {$fastaName}) {
        print OUT 
            "@" . "$fastaName\n$sequences->{$fastaName}->{'seq'}" . "+" . "$fastaName\n$sequences->{$fastaName}->{'quals'}";
    }
}

close OUT;

#############################################
sub LoadFastq {
  # read FASTQ file into a hash of hashes
  my $fastqFile=shift;
  # ensure fastq file exists else die
  open (IN, "<$fastqFile") or die ("Can not find your input file $fastqFile: No such file or directory\n");

  # set up necessary arrays
  my $names;
  my %sequences;
  my $list;

  # for as long as the file has more data:
  while (<IN>){
    chomp;
    # place every line in one of the two hashes, 
    # using the @ symbol to identify the start of each sequence
    if ($_ =~ /@/) {
	$names = $_;
	$names=~ s/@//;
	$sequences{$names} -> {'seq'} = <IN>;
	$list = <IN>;
	$sequences{$names} -> {'quals'} = <IN>;
    }
  }
  close IN;

  return \%sequences;
}

############

__END__

=head1 SYNOPSIS

Fastq_Read_Extract_byReadName.pl [Options]

Purpose: To find a list of reads based upon their read names.
The reads will be listed in the order of the queries within the
query file. Query names should not start with "@", ">", or "+".

Options:

    -help/-h     brief help message
    -input/-i    input fastq file
    -query/-q    query file [ each read name should be on a separate line ]
    -output/-o   output fastq file name

=cut
