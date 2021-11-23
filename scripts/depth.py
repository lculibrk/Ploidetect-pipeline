#!/usr/bin/env python3

import pysam
import sys
import math
import re
import collections
import argparse

def Bases_At_Pos(samfile, pos, chromname, minmapqual):
    'Return a string of the bases at that position.'
    position = 0
    coverage = 0
    bases = ""
    for pileupcolumn in samfile.pileup(reference=chromname, start=pos-1, end=pos):
        if ((pileupcolumn.pos+1) >= pos and (pileupcolumn.pos+1) <= pos):
            position = int(pileupcolumn.pos+1)
            coverage = int(pileupcolumn.n)
            for pileupread in pileupcolumn.pileups:
                if (pileupread.indel == 0 and pileupread.is_del == 0 and \
                float(pileupread.alignment.mapq) >= minmapqual):
                    bases += pileupread.alignment.seq[pileupread.query_position]
    return bases


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
    	description="Compute per-read depth",
    	formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
        parser.add_argument("-c", "--cram", required = True, help = "Path to cram")
        parser.add_argument("-o", "--output", required = True, help = "output file")
        parser.add_argument("-r", "--region", required = True, help = "region/chromosome to compute")
        parser.add_argument("-q", "--qual", required = True, help = "Minimum quality")
        parser.add_argument("-m", "--maxd", required = True, help = "max depth")
        parser.add_argument("-f", "--final", required = True, help = "Final output checkpoint file")
        args = parser.parse_args()
        ## Read cram file
        samfile = pysam.AlignmentFile(args.cram)
        chrom = args.region
        ## If output file exists, open it
        if os.path.isfile(args.output):
            with open(args.output) as f:
                out = f.readlines()
            final_line = out[-1]
            ## If last written line was truncated somehow, strip it out
            if not "\n" in final_line:
                out.pop()
                with open(args.output) as f:
                    f.writelines(out[:-1])
                final_line = out[-2]
            final_line = final_line.strip().split("\t")
            start_pos = final_line[1] + 1
        else:
            start_pos = 1

        for pileupcolumn in samfile.pileup(chrom, start_pos, min_mapping_quality = parser.qual, compute_baq = False, stepper = "nofilter", max_depth = parser.maxd):
            print(chrom + "\t" + str(pileupcolumn.pos) + "\t" + str(pileupcolumn.n))
            sys.stdout.flush()
        ## If we're done, write the final output
        with open(args.final) as f:
            f.writelines(["done\n"])
