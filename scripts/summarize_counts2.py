#!/usr/bin/env python3
"""summarize_counts.py - bins.txt
# Inputs: a stream of input ("-", stdin) from samtools depth with an "end" column added, and a bed file of intervals
# Outputs mean of stream over intervals in bins.txt
# Outputs bins to stdout
"""
import fileinput
import sys
import os

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


def getwholeline(f, byte):
    if byte > os.path.getsize(f):
        byte = os.path.getsize(f) - 2
    elif byte == os.path.getsize(f):
        byte = byte - 2
    else:
        byte = max(0, byte - 1)
    f = open(f, "rb")
    f.seek(byte, 0)
    upstream = False
    ## Search upstream for newline
    if byte > 0:
        while not upstream:
            s = f.read(1).decode('utf-8')
            if s != "\n":
                byte = byte - 1
                f.seek(byte, 0)
            else:
                byte = byte + 1
                upstream = True
                f.seek(byte, 0)
    downstream = False
    linelen = 1
    f.seek(byte, 0)
    while not downstream:
        s = f.read(linelen).decode('utf-8')
        f.seek(byte, 0)
        if s.endswith("\n"):
            downstream = True
        else:
            linelen += 1
    f.close()
    return({"line":s, "byte":byte})

def index_chromosomes(f):
    index_dict = {}
    fsize = os.path.getsize(f)
    a,b = 0,fsize
    m = round((a + b)/2)
    final_chrom = getwholeline(f, b)["line"].strip().split("\t")[0]
    chromosome_processing = True
    chrom = getwholeline(f, a)["line"].strip().split("\t")[0]
    index_dict[chrom] = a
    while chromosome_processing:
        search_chrom = getwholeline(f, a)["line"].strip().split("\t")[0]
        chrom2 = getwholeline(f, b)["line"].strip().split("\t")[0]
        searching = True
        s,e = a, b
        while searching:
            ## s = start interval, e = end interval, m = midpoint
            ## We need to search smaller ie. in (s, m)
            ## Check midpoint
            mid_line = getwholeline(f, m)
            ## Midpoint is too large, search [s,m]
            if mid_line["line"].strip().split("\t")[0] != search_chrom:
                e = m
                m = round((s + e)/2)
            ## Midpoint is too close, search [m,e]
            else:
                s = m
                m = round((s + e)/2)
            l1 = getwholeline(f, s)
            l1b = l1["byte"]
            l1 = l1["line"]
            chrom = l1.strip().split("\t")[0]
            l2 = getwholeline(f, m)["line"]
            chrom2 = l2.strip().split("\t")[0]
            ## Get line after the start to see if we're on sequential lines now
            fline = getwholeline(f, l1b + len(l1))
            if fline["line"] == l2 and chrom != chrom2:
                searching = False
        a = fline["byte"]
        m = round((a + b)/2)
        e = os.path.getsize(f)
        index_dict[fline["line"].strip().split("\t")[0]] = fline["byte"]
        if fline["line"].strip().split("\t")[0] == final_chrom:
            chromosome_processing = False
    return(index_dict)

def digit_count(n):
    k = len(str(n))
    s = k * (n + 1) - int((10**k - 1)/9)
    return(s)

def hunga_bunga(n):
    r = range(1, n + 1)
    r = [str(i) for i in r]
    r = "".join(r)
    return(len(r))

def position_calculator(index_dict, chrom, pos, dplen = 5):
    whitespace_n = 3
    offset = index_dict[chrom]
    plen = len(str(pos))
#    print(f"whitespace contrib:{whitespace_n*pos}")
#    print(f"chromosome contrib:{len(chrom)*pos}")
#    print(f"depth contrib:{dplen * pos}")
    base_line = whitespace_n + len(chrom) + dplen
#    print(digit_count(pos))
    s = offset + (base_line * pos) + digit_count(pos)
    return(s)

## Loop over bin lines
index = index_chromosomes(sys.argv[1])

#print("Position with hunga bunga method:")
#print(hunga_bunga(int(sys.argv[3])))

#print("\n")
#print("Position with mathy method:")
#print(digit_count(int(sys.argv[3])))
def bytesto(bytes, to, bsize=1024):
    """convert bytes to megabytes, etc.
       sample code:
           print('mb= ' + str(bytesto(314575262000000, 'm')))
       sample output: 
           mb= 300002347.946
    """

    a = {'k' : 1, 'm': 2, 'g' : 3, 't' : 4, 'p' : 5, 'e' : 6 }
    r = float(bytes)
    for i in range(a[to]):
        r = r / bsize

    return(r)
#print(getwholeline(sys.argv[1], position_calculator(index, "1", int(sys.argv[3]), 5)))
import re
b_limit = 10000000
halt = False
f = open(sys.argv[1])
with open(sys.argv[2]) as bins:
    for b in bins:    
        b = b.strip().split()
        chrom = b[0]
        const_len = len(chrom) + 8
        start = int(b[1])
        end = int(b[2])
        running_sum = 0
        n = 0
        linestart = position_calculator(index, chrom, start)
        lineend = position_calculator(index, chrom, end - 1)
#        print(b)
        #print(getwholeline(sys.argv[1], linestart))
        #print(getwholeline(sys.argv[1], lineend))
        f.seek(linestart, 0)
        read_start = linestart
        read_end = min([lineend, getwholeline(sys.argv[1], linestart + b_limit)["byte"]])
        while read_end != lineend:
#            print("Experimental stuff beginning")
            read_bytes = read_end - read_start
#            print(f"Supposed to end at {lineend}, but this is larger than 1Mb")
#            print(f"instead ending at byte {read_bytes} for this chunk")
            start_line = getwholeline(sys.argv[1], read_start)["line"]
            end_line = getwholeline(sys.argv[1], read_end)["line"]
#            print(f"start line of interval: \n {start_line}")
#            print(f"end line of interval: \n {end_line}")
            lines = f.readlines(read_bytes - 1)
            while lines[-1][0:len(chrom)] != chrom:
                lines.pop()
            for line in lines:
                running_sum += int(line[-6:-1])
                pos = int(line[(len(chrom)+1):-7])
                if pos < start or pos > end:
                    raise ValueError(f"Position of {pos} outside bound of ({start}, {end})")
            n += len(lines)
            read_start = read_end
            read_end = min([lineend, getwholeline(sys.argv[1], read_start + b_limit)["byte"]])
        bytes_to_read = read_end - read_start
        lines = f.readlines(bytes_to_read-1)
        while lines[-1][0:len(chrom)] != chrom:
            lines.pop()
        n += len(lines)
        for line in lines:
            running_sum += int(line[-6:-1])
            pos = int(line[(len(chrom)+1):-7])
            if pos < start or pos > end:
                raise ValueError(f"Position of {pos} outside bound of ({start}, {end})")
        #lines = [int(line[2]) for line in lines if line]
        print("\t".join(b) + "\t" + str(running_sum/n))
        sys.stdout.flush()
f.close()

sys.exit()


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
