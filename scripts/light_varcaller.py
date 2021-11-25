import pysam
import sys
import math
import re
import collections
import argparse
import os
import shutil

def Bases_At_Pos(samfile, pos, chromname, minbasequal, minmapqual):
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
                (ord(pileupread.alignment.qual[pileupread.query_position])-33) >= minbasequal and \
                float(pileupread.alignment.mapq) >= minmapqual):
                    bases += pileupread.alignment.seq[pileupread.query_position]   
    return bases

def get_reference_base(fasta, fai, chromname, pos):
        fasta.seek(0,0)
        linebyte = int(fai[chromname][3])
        linelen = int(fai[chromname][2])
        start = pos-1 
        end = pos
        startline = math.floor((start-1)/linelen)
        endline = math.floor(end/linelen)
        startbyte = start + startline
        endbyte = end + endline
        fasta.seek(int(fai[chromname][1]) + startbyte)
        seq = fasta.read(endbyte-startbyte).replace("\n", "").upper()
        return seq


if __name__ == "__main__":
        parser = argparse.ArgumentParser(
            description="Compute per-read depth",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        parser.add_argument("-n", "--ncram", required = True, help = "Path to normal cram")
        parser.add_argument("-t", "--tcram", required = True, help = "Path to tumour cram")
        parser.add_argument("-o", "--output", required = True, help = "output file")
        parser.add_argument("-r", "--regions", required = True, help = "positions file")
        parser.add_argument("-a", "--fasta", required = True, help = "Fasta file")
        parser.add_argument("-f", "--final", required = True, help = "Final output file")
        args = parser.parse_args()
        normfile = pysam.AlignmentFile(args.ncram)
        somafile = pysam.AlignmentFile(args.tcram)
        positions_file = args.regions
        fasta = open(args.fasta, "r")
        with open(args.fasta + ".fai") as f:
                fai = f.readlines()
        fai = [line.strip().split("\t") for line in fai]
        fai = {line[0]:line[1:] for line in fai}
        ## Find last written position
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
                    f.seek(pos + 1, os.SEEK_SET)
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
            final_line = final_line.strip().split("\t")
            offset = int(final_line[4])
        else:
            final_line = ""
            offset = 0
        with open(positions_file) as p:
                p.seek(offset)
                byte_count = p.tell()
                entry = p.readline().strip().split("\t")
                if len(entry) >= 2:
                        lines_to_run = True
                else:
                        lines_to_run = False
                while lines_to_run:
                        bases = Bases_At_Pos(normfile, int(entry[1]), entry[0], 20, 20).upper()
                        ref_base = get_reference_base(fasta, fai, entry[0], int(entry[1])).upper()
                        ref_count = bases.count(ref_base)
                        af = ref_count/len(bases) if bases else 1
                        if af > 0.35 and af < 0.65 and len(bases) > 10:
                                bases = Bases_At_Pos(somafile, int(entry[1]), entry[0], 20, 20).upper()
                                if len(bases) > 10:
                                        ref_count = bases.count(ref_base)
                                        af = ref_count/len(bases)
                                        print(entry[0] + "\t" + str(entry[1]) + "\t" + str(int(entry[1]) + 1) + "\t" + format(af, '.6f') + "\t" + str(byte_count))
                                        sys.stdout.flush()
                        entry = p.readline().strip().split("\t")
                        byte_count = p.tell()
                        if len(entry) < 2:
                                lines_to_run = False
                shutil.copyfile(args.output, args.final)
