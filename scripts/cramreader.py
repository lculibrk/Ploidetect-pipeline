import pysam
import sys
import math
import re
import collections



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
        samfile = pysam.AlignmentFile(sys.argv[1])
        positions_file = sys.argv[2]
        fasta = open(sys.argv[3], "r")
        with open(sys.argv[3] + ".fai") as f:
                fai = f.readlines()
        fai = [line.strip().split("\t") for line in fai]
        fai = {line[0]:line[1:] for line in fai}
        with open(positions_file) as p:
                entry = p.readline().strip().split("\t")
                if len(entry) >= 2:
                        lines_to_run = True
                else:
                        lines_to_run = False
                while lines_to_run:
                        bases = Bases_At_Pos(samfile, int(entry[1]), entry[0], 0, 0).upper()
                        ref_base = get_reference_base(fasta, fai, entry[0], int(entry[1])).upper()
                        ref_count = bases.count(ref_base)
                        af = format(ref_count/len(bases), '.6f')
                        print(entry[0] + "\t" + str(entry[1]) + "\t" + str(int(entry[1]) + 1) + "\t" + str(af) + "\t" + str(ref_base) + "\t" + str(bases))
                        sys.stdout.flush()
                        entry = p.readline().strip().split("\t")
                        if len(entry) < 2:
                                lines_to_run = False
