"""
Mappability.py
by: Tovah Markowitz
date: 12/09/15
purpose: to determine what portions of the genome map uniquely to either SK1 or S288C
input: two SAM files of the same split genome mapped to two genomes by Bowtie
output: two smaller SAM files only containing uniquely mapping portions of the genome
"""

############
#  Modules

import optparse

############
# Functions

def ParseSAM (SAMfile):
    # Purpose: to parse a SAM file into a list of header rows and
    # a dictionary of mapped reads with the keys being the read names
    # NOTE: function returns two separate outputs
    f = open (SAMfile,'r')
    header = [ ]
    reads = { }

    # since files are large I will remove empty rows as I go
    for line in f:
        # first extract header information
        if line.startswith("@"):
            header.append( line.strip('\n') )
            # if not a header, make sure it mapped
            # make read number into key, and entire row into value
        elif line.find("MD") != -1:
            reads[ line.split("\t")[0] ] = line

    return (header,reads)

def UniqueMapping (SAMfile1, SAMfile2):
    # Purpose: to identify reads that only map to one of two SAM files
    # Reads solely mapping to the first file will be called "unique"
    # Reads solely mapping to the second file will be called "error"

    # parse both SAM files
    (header1, reads1) = ParseSAM(SAMfile1)
    (header2, reads2) = ParseSAM(SAMfile2)
    
    # get all reads that were mapped to individual genome
    readName1 = set(reads1.keys())
    readName2 = set(reads2.keys())
    
    # use set difference to find reads that only mapped to one genome or the other
    # files in error didn't mapped to the original genome, but did map to the other genome
    UniqueMap = readName1.difference(readName2)
    ErrorMap = readName2.difference(readName1)
    
    # get reads that mapped to only original genome or opposite genome
    outUnique = []
    for readName in sorted(UniqueMap):
        outUnique.append( reads1.get(readName) )

    outError = []
    for readName in sorted(ErrorMap):
        outError.append( reads2.get(readName) )

    outUName = SAMfile1.split(".")[0] + "-Unique.sam"
    outEName = SAMfile2.split(".")[0] + "-Unique.sam"
    f = open(outUName,'w')
    f.write( '\n'.join(i for i in header1) )
    f.write( '\n' )
    f.write( ''.join(i for i in outUnique) )
    f.close()
    f = open(outEName,'w')
    f.write( '\n'.join(i for i in header2) )
    f.write( '\n' )
    f.write( ''.join(i for i in outError) )
    f.close()
    return (outUnique)



################################################################################
# Main

desc="""
This script is designed to compare two SAM files created by Bowtie. It identifies
reads uniquely mapped to either of the two SAM files and saves them as new SAM files.
Useful for determining which portions of the genome are unique. Output files are called 
[SAMfile1]-Unique.sam and [SAMfile2]-Unique.sam.
"""

# parse object for managing input options
parser = optparse.OptionParser(description=desc)

# essential data, defines commandline options
parser.add_option('-1', dest= 'SAMfile1', default='', help= "This input is the \
name of the 1st SAM file.")
parser.add_option('-2', dest= "SAMfile2", default='', help= "This input is the \
name of the 2nd SAM file.")

# load the inputs
(options, args) = parser.parse_args()

# read the inputs from commandline
SAMfile1 = options.SAMfile1
SAMfile2 = options.SAMfile2

UniqueMapping( SAMfile1, SAMfile2 )

#a = UniqueMapping("AH119A-053013-PM.sam","AH119A-053013-PM.sam")
