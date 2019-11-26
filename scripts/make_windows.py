import sys
import fileinput
count = 0
start = 0
infile = fileinput.input(files="-")
threshold=int(sys.argv[2])
chrom = None
for line in infile:
    # File handling and santisation
    line = line.rstrip('\n')
    line = line.split("\t")
    if not chrom:
        chrom = line[0]
    # Check if new chromosome
    if(chrom != line[0]):
        # Output last chromosome, previous end (start), end of chromosome (lastknownend), count
        out = [chrom, start, lastknownend, count]
        # Convert out to strings and join them for printing, print
        out = [ str(x) for x in out ]
        out = "\t".join(out)
        print(out)
        # Get things ready for next loop
        chrom = line[0]
        start = int(line[1])
        count = 0
    if(count < threshold):
        # Add to count
        count += int(line[3])
        # Record end position of the count (in case new chromosome is skipped to)
        lastknownend = int(line[2])
    if(count >= threshold):
        stop = int(line[2])
        chrom = line[0]
        out = [chrom, start, stop, count]
        out = [ str(x) for x in out ]
        out = "\t".join(out)
        start = int(line[2])
        count = 0
        print(out)
    sys.stdout.flush()
