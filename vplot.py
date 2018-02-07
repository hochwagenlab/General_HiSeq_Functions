# vplot.py
# made by Tovah Markowitz
# uses python3
# 1/22/18

#########################
# MODULES
import matplotlib.pyplot as plt
import numpy
import optparse
from datetime import datetime
import pysam # requires 2.7.12 on HPC

#########################
# FUNCTIONS

def import_sam_pysam(sam_file):
    """Import sam file using pysam and save only query name, chromosome
    name, start position, and fragment length as a list of lists.
    """
    print("Reading in SAM file.")
    print(datetime.now())
    sam_in = pysam.AlignmentFile(sam_file,"r")
    sam_in2 = [row for row in sam_in.fetch(until_eof=True) ]    
    sam_map = [[row.query_name,
                row.reference_name, # chromosome name
                (row.reference_start + 1), # 0-based
                row.template_length ] # fragment length
                   for row in sam_in2 if (row.flag != 77) & (row.flag != 141)]
    print("Finished.")
    print(datetime.now())
    return sam_map

def parse_pair_pysam(sam_list):
    """Parse paired-end sequencing data to extract chromosomes, 
    fragment midpoint, and fragment length. The input should be the output
    of import_sam_pysam. The output is a dictionary where the keys are 
    chromosomes and the values are lists of midpoint and fragment length 
    pairs. This function only works when reads can only be mapped once.
    """
    sam_list.sort()
    chrs = set(row[1] for row in sam_list)
    sam_dict = {chrom: [ ] for chrom in chrs}
    for i in range(len(sam_list) - 1):
        if (i % 5000000) == 0:
            print("Parsed " + str(i) + " lines.")
            print(datetime.now())
        if " " in sam_list[i][0]:
            print("Caution: read was not processed. " + sam_list[i][0]) 
        elif sam_list[i][0] == sam_list[i + 1][0]:
            if sam_list[i][1] == sam_list[i + 1][1]:
                    chrom = sam_list[i][1]
                    frag_len = abs(int(sam_list[i][3]))
                    start = min(int(sam_list[i][2]), int(sam_list[i + 1][2]))
                    end = start + frag_len
                    mid = int((start + end) / 2)
                    sam_dict[chrom].append([mid, frag_len])
            else:
                print("Caution: paired reads mapped to different chromosomes. "
                          + sam_list[i][0])
    return sam_dict
        
def import_sam(sam_file):
    """Read in sam file of sequencing data. The output is a list of lists 
    containing data on all mapping reads.
    """
    print("Reading in SAM file.")
    print(datetime.now())
    f = open(sam_file, 'r')
    sam_in = f.readlines()
    f.close
    sam = [row.strip().split('\t') for row in sam_in
               if not row.startswith("@")]        
    # removing rows with FLAG indicating no mapping of paired reads
    sam_map = [row for row in sam if (row[1] != '77') & (row[1] != '141')]
    print("Finished.")
    print(datetime.now())
    return sam_map

def parse_pair_sam(sam_list):
    """Parse paired-end sequencing data to extract chromosomes, fragment 
    midpoint, and fragment length. The input should be the output of 
    import_sam. The output is a dictionary where the keys are chromosomes 
    and the values are lists of midpoint and fragment length pairs. This 
    function only works when reads can only be mapped once.
    """
    sam_list.sort()
    chrs = set(row[2] for row in sam_list)
    sam_dict = {chrom: [ ] for chrom in chrs}
    for i in range(len(sam_list) - 1):
        if (i % 5000000) == 0:
            print("Parsed " + str(i) + " lines.")
            print(datetime.now())
        if " " in sam_list[i][0]:
            print("Caution: read was not processed. " + sam_list[i][0]) 
        elif sam_list[i][0] == sam_list[i + 1][0]:
            if sam_list[i][2] == sam_list[i + 1][2]:
                if sam_list[i][3] == sam_list[i + 1][7]:
                    chrom = sam_list[i][2]
                    frag_len = abs(int(sam_list[i][8]))
                    start = min(int(sam_list[i][3]), int(sam_list[i][7]))
                    end = start + frag_len
                    mid = int((start + end) / 2)
                    sam_dict[chrom].append([mid, frag_len])
                else:
                    print("Caution: pair positions don't match. "
                              + sam_list[i][0])
            else:
                print("Caution: paired reads mapped to different chromosomes. "
                          + sam_list[i][0])
    return sam_dict

def extract_lengths(sam_dict):
    """Extract fragment lengths of all mapped reads from dictionary.
    """
    lengths = [ ]
    for key in sam_dict.keys():
        tmp = [row[1] for row in sam_dict[key]]
        lengths.extend(tmp)
    lengths.sort()
    return lengths

def plot_hist(data, out_hist, title="", xlab=""):
    """Plot a list of data as a histogram.
    """
    if "." in out_hist:
        outfile_name = out_hist
    else:
        outfile_name = out_hist + ".pdf"
#    plt.hist(data, bins=100, density=True)
    plt.hist(data, bins=100, normed=1)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel("Density")
    plt.savefig(outfile_name)
    plt.close("all")

def import_bed(bedfile):
    """Read in bed file and extract chromosome and midpoints.
    Orientation ignored.
    """
    f = open(bedfile, 'r')
    bed_in = f.readlines()
    f.close
    bed = [row.strip().split('\t') for row in bed_in]
    chrs = set(row[0] for row in bed)
    midpoints = {chrom: [ ] for chrom in chrs}
    for row in bed:
        midpoint = int((int(row[1]) + int(row[2]))/2) 
        midpoints[row[0]].append(midpoint)
    for key in midpoints.keys():
        midpoints[key].sort()
    return midpoints

def find_nearby(sam_dict, midpoints, max_range):
    """Identify all reads that are within a certain distance of a feature
    of interest. Assumes both inputs [sam_dict and midpoints] are 
    dictionaries with keys being chromosome names. Max_range is the maximum 
    distance from the midpoints that will be used on the plot.
    """
    nearby = [ ]
    for key in midpoints.keys():
        print("Starting to process " + key + " for vplot")
        print(datetime.now())
        for midpoint in midpoints[key]:
                if sam_dict.get(key):
                        tmp=[[read[0] - midpoint, read[1]] for read
                                     in sam_dict[key]
                                     if abs(midpoint - read[0]) < max_range]
                        nearby.extend(tmp)
    return nearby

def group_data(nearby, dist_window_size, fraglen_window_size,
                   max_range,max_frag_length):
    """Group data from find_nearby into windows based upon both fragment
    length and distance from feature midpoints.
    """
    maxfrag = max([row[1] for row in nearby])
    if maxfrag > max_frag_length:
        maxfrag = max_frag_length
    minfrag = min([row[1] for row in nearby])
    frag_cut = [a for a in range(minfrag, maxfrag, fraglen_window_size)]
    dist_cut = [a for a in range(
        int(max_range / dist_window_size) * -dist_window_size,
        max_range + dist_window_size, dist_window_size)]
    grouped_data = [[] for _ in range(len(frag_cut)-1)]
    for i in range(len(frag_cut)-1):
        within_frag = [data for data in nearby if
                           frag_cut[i] <= data[1] < frag_cut[i+1]] 
        for j in range(len(dist_cut)-1):
            within_dist = [data for data in within_frag if
                               dist_cut[j] <= data[0] < dist_cut[j+1]]
            grouped_data[i].append(len(within_dist))
    return(grouped_data, frag_cut, dist_cut)

def vplot(grouped_data, frag_cut_points, dist_cut_points, out_vplot, xlab="",
              title=""):
    """Plot group_data as a heatmap with axes identifying distances.
    """
    if "." in out_vplot:
        outfile_name = out_vplot
    else:
        outfile_name = out_vplot + ".pdf"
    plt.imshow(grouped_data, cmap=plt.cm.bone_r, origin='lower', aspect='auto')
    plt.colorbar()
    min_dist = dist_cut_points[0]
    max_dist = dist_cut_points[-1]
    byX = len(grouped_data[0]) / 4
    plt.xticks(numpy.arange(-0.5, len(grouped_data[0]), byX),
                   [min_dist, min_dist / 2, 0, max_dist / 2, max_dist])
    min_frag = frag_cut_points[0]
    max_frag = frag_cut_points[-1]
    frag_dist = (max_frag - min_frag) / 4
    byY = len(grouped_data) / 4
    plt.yticks(numpy.arange(-0.5, len(grouped_data), byY),
                   [min_frag, min_frag + frag_dist, min_frag + 2*frag_dist,
                    max_frag - frag_dist, max_frag])
    plt.title(title)
    if xlab:
        plt.xlabel(xlab)
    else:
        plt.xlabel("Distance from feature (bp)")
    plt.ylabel("Fragment size (bp)")
    plt.savefig(outfile_name)
    plt.close("all")

#########################
# MAIN

    
desc="""
This script is designed to replicate that of the v-plot used by the Henikoff
lab. See PMID: 22025700 for details. Specifically, the plot is a 2-D heatmap
where the x-axis is the distance from a specified feature (such as a motif) and
the y-axis is the fragment size, and intensity indicates the number of reads
within a given bin. In the words of the Henikoff lab: the length of each 
fragment was plotted as a function of the distance from the fragment midpoint 
to the center of the site for each annotated feature. Heatmaps were generated 
by quantifying the number of reads at each relative distance and length 
coordinate.
"""

# parse object for managing input options
parser = optparse.OptionParser(description=desc)

# essential data, defines commandline options
parser.add_option('-b', dest='bedfile', default='', help="The bed file \
fragments should be compared to. Optional: only required for vplot.")
parser.add_option('--hist', dest='outhist', default='', help="The name of \
the histogram output file. Optional: Necessary only for making histogram.")
parser.add_option('--maxdist', dest="maxdistance", default=1000, help="Maximum \
distance away from feature of interest to be analyzed. Determines x-axis of \
vplot. Default: 1000 bp")
parser.add_option('--maxlen', dest="maxlength", default=500, help="Maximum \
sequenced fragment length analyzed. Default: 500 bp")
parser.add_option('-s', dest='samfile', default ='', help="The SAM file to \
be analyzed.")
parser.add_option('-t', dest='title', default='', help="The title for the \
figures. Optional.")
parser.add_option('-v', dest="outvplot", default='', help="The name of the \
vplot outplot file. Optional: only required if making a vplot.")
parser.add_option('--windist', dest='windowdist', default=25, help="Window \
size for x-axis (distance from feature) of vplot. Default: 25 bp")
parser.add_option('--winlen', dest='windowlen', default=5, help="Window \
size for y-axis (fragment size) of vplot. Default: 5 bp")
parser.add_option('-x', dest='xlabvplot', default='', help="The label for the \
x-axis on the vplot. For example, 'distance from hotspots'. Optional.")

# load the inputs
(options,args) = parser.parse_args()

# reads the inputs from commandline
bedfile = options.bedfile
maxdistance = options.maxdistance
maxlength = options.maxlength
outhist = options.outhist
outvplot = options.outvplot
samfile = options.samfile
title = options.title
windowdist = options.windowdist
windowlen = options.windowlen
xlabvplot = options.xlabvplot

if __name__ == "__main__":
    sam_list = import_sam_pysam(samfile)
    sam_dict = parse_pair_pysam(sam_list)
#    sam_list = import_sam(samfile)
#    sam_dict = parse_pair_sam(sam_list)
    if (outlist is not "") | (outhist is not ""):
        lengths = extract_lengths(sam_dict)
    if outhist:
        plot_hist(lengths, out_hist=outhist, title=title,
              xlab="sequenced fragment length (bp)")
    if outvplot:
        if not bedfile:
            print("No bed file indicated for vplot analysis.")
        else:
            midpoints = import_bed(bedfile)
            nearby = find_nearby(sam_dict, midpoints, int(maxdistance))
            (grouped, frag_cut_points, dist_cut_points) = group_data(nearby,
                  int(windowdist), int(windowlen), int(maxdistance),
                  int(maxlength))
            vplot(grouped, frag_cut_points, dist_cut_points, outvplot,
                       xlabvplot, title)
            



