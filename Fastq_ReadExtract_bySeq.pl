#!/usr/bin/perl                                                                                          

# Name: Fastq_ReadExtract_bySeq.pl
# Author: Tovah Markowitz  
# Date: 11/7/2013                                                                                                   
# Description: extract specific queries from a fastq file based upon sequence                                      
# Input: name of fastq file  
# Output: queries in fastq format                                                                                                                                  
use strict;
use Data::Dumper;
use Getopt::Long;                                                                             
use Pod::Usage;

my $help=0;
my ($fastqFile, $querySeq, $queryPos, $outputFile);

# set commandline options.                                        
GetOptions ("input=s" => \$fastqFile,
            "qSeq=s" => \$querySeq,
	    "qPos=i" => \$queryPos,
            "output=s" => \$outputFile,
	    "help" => \$help,) or pod2usage(2);
pod2usage(1) if $help;

#############################################

my $sequences=LoadFastq($fastqFile);

# create regular expression pattern based upon query position and query sequence
my $regexp;
if ($queryPos) {
  $regexp=join( "", "^[ACTGN]{", $queryPos-1, "}", $querySeq );
} else {
  $regexp=$querySeq;
}

open (OUT, ">$outputFile") or die ("Can not open $outputFile for writing.\n");
# for each sequence check if querySeq exists at the noted position
foreach my $name (sort keys ( %{$sequences} )) {
  if ($sequences -> {$name} -> {'seq'} =~ m/$regexp/) {
     my $at = $name;
     my $plus = $name;
     $plus =~ s/>/+/;
     $at =~ s/>/@/;
     print OUT
       "$at\n$sequences->{$name}->{'seq'}$plus\n$sequences->{$name}->{'qual'}";
  }
}
close OUT;

#############################################
sub LoadFastq {
  # read FASTQ file into a hash of hashes
  my $fastqFile = shift;
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
	$names =~ s/@/>/;
	$sequences{$names} -> {'seq'} = <IN>;
	$list = <IN>;
	$sequences{$names} -> {'quals'} = <IN>;
    }
  }
  close IN;

  return \%sequences;
}

###########################################

__END__

=head1 SYNOPSIS

Fastq_ReadExtract_bySeq.pl [Options]

Purpose: To find reads from a FASTQ file with a specific query sequence.
The function will only find reads with an exact match to the entire query. 
Calling qPos will search for the query string at a specific position within
the read.

Options:

     -help/-h        brief help message
     -input/-i       input fastq file
     -qSeq/-qS       query sequence
     -qPos/-qP       (optional) query position
     -output/-o      output fastq file name

=cut
