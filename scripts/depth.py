#!/usr/bin/env python3

import pysam
import sys
import math
import re
import collections
import argparse
import os
import shutil

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
        fsize = os.stat(args.output).st_size
        if fsize > 0:
            with open(args.output, "rb+") as f:
                ## strip final line in case it's truncated
                f.seek(0, os.SEEK_END)
                pos = f.tell() - 1
                while pos > 0 and f.read(1).decode() != "\n":
                    pos -= 1
                    f.seek(pos, os.SEEK_SET)
                if pos > 0:
                    f.seek(pos, os.SEEK_SET)
                    f.truncate()
            with open(args.output, "rb") as f:
                f.seek(0, os.SEEK_END)
                pos = f.tell() - 1
                b_read = 1
                while pos > 0 and f.read(1).decode() != "\n":
                    pos -= 1
                    f.seek(pos, os.SEEK_SET)
                    b_read += 1
                pos += 1
                final_line = f.read(b_read).decode()
            start_pos = int(final_line.split()[1]) + 1
        else:
            start_pos = 1
        caught_up = False
        for pileupcolumn in samfile.pileup(chrom, max(0, start_pos-150), min_base_quality = 0,  min_mapping_quality = int(args.qual), max_depth = int(args.maxd), compute_baq = False, ignore_orphans = False):
            if not caught_up:
                if pileupcolumn.reference_pos < start_pos:
                    continue
                else:
                    caught_up = True
            print(chrom + "\t" + str(pileupcolumn.reference_pos) + "\t" + str(pileupcolumn.n))
            sys.stdout.flush()
        ## If we're done, write the final output
        shutil.copyfile(args.output, args.final)
