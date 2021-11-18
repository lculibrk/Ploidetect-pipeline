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
infile = fileinput.input(files=sys.argv[1])
# Read bins file
with open(sys.argv[2]) as f:
    bins = f.readlines()
bins = [bin.strip().split("\t") for bin in bins]
bin = 0
chrom = bins[bin][0]
start = bins[bin][1]
end = bins[bin][2]
# Loop over input lines
for line in infile:
    # File handling and santisation
    line = line.rstrip("\n")
    line = line.split("\t")
    count = count + int(line[COUNT])
    nbins = nbins + 1
    ## Check if end of interval
    if line[START] == end:
        m = count/float(nbins)
        outline = str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(m)
        print(outline)
        sys.stdout.flush()
        if bin < len(bins) - 1:
            bin = bin + 1
            chrom = bins[bin][0]
            start = bins[bin][1]
            end = bins[bin][2]
            nbins = 0
            count = 0
    # Ensures un-buffered output to stdout
    sys.stdout.flush()
