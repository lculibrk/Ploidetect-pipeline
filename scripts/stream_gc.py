#!/usr/bin/env python3
"""stream_gc.py in.fa bins.txt
# Inputs: a fasta input (in.fa), and a bed file of intervals
# Outputs mean of GC content over intervals in bins.txt
# Outputs GC contents to stdout
"""
import sys
import math
import re
from decimal import Decimal, ROUND_HALF_UP

def normal_round(n):
    if n - math.floor(n) < 0.5:
        return math.floor(n * 1000000)/1000000
    return math.ceil(n * 1000000)/1000000

binfile = open(sys.argv[2], "r")
fasta = open(sys.argv[1], "r")
with open(sys.argv[1] + ".fai") as f:
	fai = f.readlines()

fai = [line.strip().split("\t") for line in fai]
for line in range(len(fai)):
    fai[line][0] = re.sub("chr", "", fai[line][0])
fai = {line[0]:line[1:] for line in fai}


## move pointer to start
fasta.seek(0,0)
for line in binfile:
	line = line.strip().split("\t")
	chrom = line[0]
	start = int(line[1])
	end = int(line[2])
	linebyte = int(fai[chrom][3])
	linelen = int(fai[chrom][2])
	## find lines spanning interval, N of newline bytes
	startline = math.floor(start/linelen)
	endline = math.floor(end/linelen)
	startbyte = start + startline
	endbyte = end + endline
	fasta.seek(int(fai[chrom][1]) + startbyte)
	seq = fasta.read(endbyte-startbyte).replace("\n", "").upper()
	seqlen = len(seq)
	g = seq.count("G")
	c = seq.count("C")
	## Not technically proper, but to maintain consistency with bedtools nuc we round traditionally instead of banker's rounding
	gc = format(round((g+c)/seqlen, 6), '.6f')
	print(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(gc))
	sys.stdout.flush()

