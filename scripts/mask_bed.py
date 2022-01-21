#!/usr/bin/env python3
"""mask_bed.py data.bed mask.bed
# Inputs: a data file (data.bed) and a mask file (mask.bed)
# Outputs data file with regions overlapping mask.bed set to zero
# Outputs data to stdout
"""

import fileinput
import sys
import os
import pathlib

sys.path.append(pathlib.Path(__file__).parent.resolve())

import summarize_counts2
print(__name__)
if __name__ == "__main__":
	index = summarize_counts2.index_chromosomes(sys.argv[1])
	fsize = os.stat(sys.argv[1]).st_size
	with open(sys.argv[2], "r") as f:
		print("Opened file")
		n_reps = 0
		for rep in f:
			n_reps = n_reps + 1
			rep = rep.strip().split("\t")
			startpos = int(rep[1])
			endpos = int(rep[2])
			chrom = rep[0]
			## Get the end byte
			end_byte = summarize_counts2.position_calculator(index, rep[0], endpos, 5)
			if end_byte > fsize:
				start_byte = summarize_counts2.position_calculator(index, rep[0], startpos, 5)
				if start_byte > fsize:
					print("Exiting prematurely")
					print(f"Reps completed: {n_reps}")
					sys.exit()
			while end_byte > fsize:
				endpos = endpos + 1
				end_byte = summarize_counts2.position_calculator(index, rep[0], endpos, 5)
			## Open file to edit in-place
			bed = open(sys.argv[1], "r+b")
			## Get byte position of first line to replace
			bed.seek(summarize_counts2.position_calculator(index, rep[0], int(startpos) - 1, 5))
			for i in range(startpos, endpos+1):
				## Assemble replacement line
				to_write = chrom + "\t" + str(i) + "\t00000\n" 
				## Read original line
				line = bed.read(len(to_write)).decode('utf-8')
				## Rewind the pointer move from reading
				bed.seek(-len(line), 1)
				## Split so we can verify the positions are what we expect
				t_rep = to_write.split("\t")
				t_lne = line.split("\t")
				if t_rep[1] != t_lne[1]:
					if t_rep[0] != t_lne[0] or not line.endswith("\n"):
						continue
					else:
						print(t_rep)
						print(t_lne)
						sys.exit()
				## If passed, overwrite the relevant lines
				bed.write(to_write.encode('utf-8'))
			bed.close()