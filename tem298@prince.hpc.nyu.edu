"""
ChromosomeWalk.py
by: Tovah Markowitz
date: 5/11/16
based upon: ChromosomeWalk.py
purpose: to determine the sequence of a genome directly downstream of a known starting point using fastq data
input: start sequence, length to walk, and fastq file
output: raw sequence of length chosen beginning with start sequence
"""

############
# Modules

#from Bio import SeqIO
import re
from datetime import datetime
import optparse

############
# Functions

def Walk (fastqFile, start, output):
    # to walk along a chromosome, starting with "start" sequence
    # and continuing up to 500bp away using sequencing read information

    print fastqFile
    # read in all sequences from FASTQ file
    reads = seq_from_fastq(fastqFile)
    readLen = len(reads[0])
    tmp = [ ]

    # input sequence (start) defines length of overlap
    overlap = len(start)
    print start

    # get sequence list to include reverse and complement of all reads
    reads3 = [get_reverse_complement(ind_read) for ind_read in reads]
    reads = reads + reads3
    out = [ ]
    max_len = overlap
    max_len1 = overlap+1
    i = 0

    # ways this function prevents neverending loops:
    # 1) maximum output fragment length is kept to 500bp
    # 2) if fragment length stops growing, halt loop
    # 3) if number of sequences to analyze becomes too great
    while ((max_len < 500) & (max_len < max_len1) & (len(tmp) <= 25) ):
        if i == 0:
            tmp = [start]
            print "Round:" + str(i)
            print datetime.now()
            print "Max length: " + str(max_len)
            print "Going to analyse " + str(len(tmp)) + " sequences"
        else:
            # tmp == sequences to continue analyzing/appending to
            max_len = max_len1
            print "Round: " + str(i)
            print datetime.now()
            print "Max length: " + str(max_len)
            print "Going to analyze " + str(len(tmp)) + " sequences"
        tmp3 = [ ]
        for seq in tmp:
            # tmp2 == for searching on a single sequence
            tmp2 = [ ]
            for ind_read in reads:

                # for each sequence in the searchable list tmp take last bases of length overlap as new search string
                a = re.search(seq[-overlap:],ind_read)
                if bool(a):
                    # if sequence and read have any sequence upstream of searched sequence
                    # make sure overlapping sequence is identical
                    if overlap < len(ind_read[:a.end()]):
                        if len(seq) > len(ind_read[:a.end()]):
                            if seq.endswith(ind_read[:a.end()]):
                                # double check that length of addition is within reason
                                if len(ind_read[a.end():]) <= (readLen - overlap):
                                    tmp2.append(seq + ind_read[a.end():] )
                                else: print ind_read + " has extension value that is too high"
                        else:
                            if ind_read[:a.end()].endswith(seq):
                                # double check that length of addition is within reason
                                if len(ind_read[a.end():]) <= (readLen - overlap):
                                    tmp2.append(seq + ind_read[a.end():] )
                                else: print ind_read + " has extension value that is too high"
                    else:
                        # double check that length of addition is within reason
                        if len(ind_read[a.end():]) == (readLen - overlap):
                            tmp2.append(seq + ind_read[a.end():])
                        else: print ind_read + "does not fit criteria"

                        
            # if none of the reads have the search sequence send sequence to out
            # otherwise send to tmp3
            if bool(tmp2) == False:
                out.append(seq)
            else:
                tmp3[len(tmp3):len(tmp3)] = tmp2
        
        # tmp4 == shorter fragments that are already encompassed into a longer fragment
        tmp3 = set(tmp3)
        tmp3 = list(tmp3)
        tmp4 = [ ]
        for j in range(len(tmp3)):
            for k in range(len(tmp3)):
                if len(tmp3[j]) > len(tmp3[k]): 
                    # if smaller fragment is found within larger fragment, we don't want it as an output
                    if tmp3[j].find(tmp3[k]) != -1:
                        tmp4.append(tmp3[k])
            if len(tmp3[j]) > max_len:
                   max_len1 = len(tmp3[j])
        # remove duplicates, and prepare for later steps
        tmp4 = set(tmp4)
        tmp3 = set(tmp3)
        tmp = tmp3.difference(tmp4)
        i += 1
        
    out[len(out):len(out)]= tmp3.difference(tmp4)
    j = 1
    f = open(output,'w')
    for k in out:
        f.write('>' + str(j) + '\n' + k + '\n')
        j += 1
    f.close()

def seq_from_fastq(fastqFile):
    # get all sequences from the fastq file
    f = open (fastqFile,'r')
    reads = [ ]
    for i,line in enumerate(f):
        if i % 4 == 1:
            reads.append( line.strip('\n') )
    f.close()
    return reads

def get_complement(sequence):
    conversion_map = {'A':'T' , 'G':'C' , 'C':'G' , 'T':'A', 'N':'N'}
    return ''.join( [conversion_map[i] for i in sequence] )

def get_reverse_complement(sequence):
    return get_complement( sequence )[::-1]


##################
# MAIN

desc="""
This script is designed to use FASTQ sequence data to walk along a short portion of the genome.
This function does not take into account read quality or the number of times a read is found.
It works by only searching for perfect matches of fragments (all search terms are equal length
to original start sequence). It does require: 2X the size of the FASTQ file as memory and a long
run time (4h average). To minimize neverending loops, I have incorporated the following
terminators to the complete the function: a) no additional reads are found to continue a fragment,
b) longest output fragment reaches beyond 500bp, and c) many potential sequences are identified (over 25).
"""

# parse object for managing input options
parser = optparse.OptionParser(description=desc)

# essential data, defines commandline options
parser.add_option ('-f', dest = 'fastqFile', default = '', help = "This input \
is the name of the fastq file to be used.")
parser.add_option ('-s', dest = 'start', default = '', help = "This input \
is the starting sequence to walk from. Must be all caps.")
parser.add_option ('-o', dest = 'output', default = '', help = "This input \
is the name of the newly created output file.")

# load the inputs
(options,args) = parser.parse_args()

# reads the inputs from commandline
fastqFile = options.fastqFile
start = options.start
output = options.output

Walk(fastqFile, start, output)
