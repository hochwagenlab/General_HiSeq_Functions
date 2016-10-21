"""
Name: Bedgraph2VariableStepWiggle.py
Created by: Tovah Markowitz
Date: 4/11/16
Update: 4/20/16
Change: No longer creates all chromosome file, also exits if output folder already exists

Note: currently only works with span=1
Key detail 1: Bedgraphs use a 0-based coordinate system, while wiggle files use a 1-based coordinate system.
Key detail 2: The end position of a bedgraph is not inclusive.
"""

##################################################################################
# Modules
from collections import defaultdict
import os
import sys
import optparse
from datetime import datetime

##################################################################################
# Functions

def read_bedgraph(bedgraph):
    # Purpose: to read in a bedgraph file and create a sorted dict composed of lists of lists
    # step 1: read in bedgraph
    f = open(bedgraph, 'r')
    bedG = f.readlines()
    f.close()

    # step 2: organize bedgraph into a dict (chromosomes)
    # of lists (each row of bedgraph) of lists (start, end, score)
    for i in range( len(bedG) ):
        bedG[i] = bedG[i].strip().split('\t')

    bedD = defaultdict(list)
    for i in bedG:
        bedD[i[0]].append( map(float, i[1:]) )

    # step 3: for each chromosome do numeric sort by start
    for key in bedD.keys():
        bedD[key].sort()
    
    return bedD

def create_variable_wiggle(bedgraph, bedD):
    # Purpose: to convert a bedgraph to a variable-step wiggle file
    # step 1: get information from file name
    filename = bedgraph.split(".")[0]
    location = os.getcwd()
    if os.path.isdir(filename):
        exit()    
    else:
        os.mkdir(filename)
        os.chdir(filename)

    # if bedgraph file includes SPMR, it means it is from MACS2 with SPMR analysis
    if filename.find("SPMR") != -1:
        dataSource = "Extended tag pileup from MACS v2.1.0 with SPMR normalization for every 1 bp"
    else:
        dataSource = filename
    # create all chromosome file if it doesn't exist and empty if it does
#    g = open( filename + "_all.wig", 'w')
#    g.close()
        
    # step 2: determine chromosome naming system
    # chromosome list for SK1 and SacCer3 in chromosome number order
    SK1 = ( "chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08",
            "chr09", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16" )
    SacCer3 = ( "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
                "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI" )
    # use first key as method to determine which genome data was mapped to
    if bedD.keys()[0] in SK1:
        chromSet = SK1
    elif bedD.keys()[0] in SacCer3:
        chromSet = SacCer3
    else:
        print "Do not recognize chromosome names.\n"
        sys.exit()

    # step 3: for each chromosome make header, do calculations, and write to file
    for chrName in chromSet:
        print datetime.now().ctime()
        # create header
        header = ( "track type=wiggle_0 name=" + filename + "_" + chrName + " description=" + dataSource +
                   "\nvariableStep chrom=" + chrName + " span=1\n" )

        out = [ ]
        # for each row in bedgraph for individual chromosome
        for row in range( len( bedD[chrName] ) ):
            # ensure that end of previous row is not bigger than start of current row or trim start position of row
            if row != 0:
                if bedD[chrName][row][0] < bedD[chrName][row-1][1]:
                    print ( "Warning: Overlaps begin at " + chrName + ":" + str( int( bedD[chrName][row][0] ) ) +
                            "-" + str( int( bedD[chrName][row][1] ) ) + " ... Trimming row\n" )
                    bedD[chrName][row][0] = bedD[chrName][row-1][1]
            # skip rows with zero as score
            if bedD[chrName][row][2] != 0:
            # expand positions and convert to wiggle numbering system and add score for each position
                positions = range( int( bedD[chrName][row][0] ) + 1, int( bedD[chrName][row][1] ) + 1 )
                tmp = "\n".join( [ str(position) + "\t" + str( bedD[chrName][row][2] ) for position in positions ] )
                out.append(tmp)
                
        # to only print out files for chromosomes with information
        if ( len(out) != 0 ):
            # write wig file for individual chromosome
            f = open( filename + "_" + chrName + ".wig", 'w')
            f.write( header )
            f.write( '\n'.join(out) )
            f.close()
            
        # write to wig file for all chromosomes
 #       g = open( filename + "_all.wig", 'a')
 #       g.write( header )
 #       g.write( '\n'.join(out) )
 #       g.write( '\n' )
 #       g.close()

    os.chdir(location)

################################################################################
# Main

desc="""
A script to convert from a bedgraph into a variableStep wiggle file.
It is NOT designed to work with bedgraphs with overlapping fragments.
Creates wiggle files for each individual chromosome with information.
Uses input file name to determine output file names.
Note: this function can currently only handle span=1.
"""

# parse object for managing input options
#parser = optparse.OptionParser()
parser = optparse.OptionParser(description=desc)

# essental data, defines commandline options
parser.add_option('-b', dest= "bedgraph", default= '', help= "This is the name \
 of the input bedgraph file.")

# load the inputs
(options,args) = parser.parse_args()

# reads the inputs from commandline
bedgraph = options.bedgraph

bedgraphDict = read_bedgraph( bedgraph )
a = create_variable_wiggle( bedgraph, bedgraphDict )
