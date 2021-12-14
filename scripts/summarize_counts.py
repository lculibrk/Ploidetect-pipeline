#!/usr/bin/env python3
"""summarize_counts.py - bins.txt
# Inputs: a stream of input ("-", stdin) from samtools depth with an "end" column added, and a bed file of intervals
# Outputs mean of stream over intervals in bins.txt
# Outputs bins to stdout
"""
import fileinput
import sys

[CHROM, START, COUNT] = range(3)

# Initialize count and start variables
count = 0
nbins = 0
start = stop = lastknownend = 0
# Read infile
infile = sys.argv[1]
# Read bins file
binfile = open(sys.argv[2])
bin = binfile.readline().strip().split("\t")
chrom = bin[0]
start = int(bin[1])
end = int(bin[2])
# Loop over input lines
with open(infile) as f:
    for line in f:
        if len(line) >= 3:
            # File handling and santisation
            line = line.rstrip("\n")
            line = line.split("\t")
            ## Check if chromosome exceeded                                                                                                                                                                                                    
            if line[CHROM] != chrom or int(line[START]) > end:
                if nbins > 0:
                    m = count/float(nbins)
                    outline = str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(m)
                    print(outline)
                    sys.stdout.flush()
                    count = int(line[COUNT])
                    nbins = 1
                bin = binfile.readline().strip().split("\t")
                if bin:
                    chrom = bin[0]
                    start = int(bin[1])
                    end = int(bin[2])
                continue
            count = count + int(line[COUNT])
            nbins = nbins + 1
            ## Check if end of interval
            if int(line[START]) == end:
                m = count/float(nbins)
                outline = str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(m)
                print(outline)
                sys.stdout.flush()
                bin = binfile.readline().strip().split("\t")
                if bin:
                    chrom = bin[0]
                    start = int(bin[1])
                    end = int(bin[2])
                    nbins = 0
                    count = 0
binfile.close()
