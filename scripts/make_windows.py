#!/usr/bin/env python3
#
# make_windows.py
# Inputs: a stream of input ("-", stdin) from samtools and a threshold (integer)
# Sums stream of input depth until cumulative coverage exceeds threshold, and makes a new bin when that occurs
# Outputs bins to stdout
#
#
# Import modules
import fileinput
import sys

# Initialize count and start variables
count = 0
start = 0
# Read infile
infile = fileinput.input(files="-")
# Get threshold
threshold = int(sys.argv[2])
# Initialize chromosome as None
chrom = None
# Loop over input lines
for line in infile:
    # File handling and santisation
    line = line.rstrip("\n")
    line = line.split("\t")
    # If first iteration, set chrom to the first chromosome value
    if not chrom:
        chrom = line[0]
    # Check if new chromosome
    if chrom != line[0]:
        # Output last chromosome, previous end (start), end of chromosome (lastknownend), count
        out = [chrom, start, lastknownend, count]
        # Convert out to strings and join them with tabs
        out = [str(x) for x in out]
        out = "\t".join(out)
        # Print to stdout
        print(out)
        # Set new chromosome
        chrom = line[0]
        # Set new start
        start = int(line[1])
        # STOP THE COUNT!
        count = 0
    # Check if all is normal and we continue with creating the bin
    if count < threshold:
        # Add to count
        count += int(line[3])
        # Record end position of the count (in case new chromosome is skipped to)
        lastknownend = int(line[2])
    # If count exceeds threshold
    if count >= threshold:
        # Record end position
        stop = int(line[2])
        # Record chromosome
        chrom = line[0]
        # Output chromosome, start of bin, end of bin, count
        out = [chrom, start, stop, count]
        # Convert out to strings and join them with tabs
        out = [str(x) for x in out]
        out = "\t".join(out)
        # Set new start
        start = int(line[2])
        # STOP THE COUNT!
        count = 0
        # Print to stdout
        print(out)
    # Ensures un-buffered output to stdout
    sys.stdout.flush()
