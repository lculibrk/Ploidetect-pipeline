import sys


## 1st arg is baf file
## 2nd arg is bin file

## Start by opening the windows file and loading the first bin

bins = open(sys.argv[2], "r")
bin = bins.readline().strip().split("\t")
[chrom, start, end, na] = bin

baflist = []
with open(sys.argv[1]) as data:
    for line in data:
        line = line.strip().split("\t")
        ## Reached the end of a bin
        while chrom != line[0] or int(line[2]) > int(end):
            ## Print out aggregated BAF
            if baflist:
                print("\t".join([chrom, start, end, ";".join(baflist)]))
            else:
                print("\t".join([chrom, start, end, "NA"]))
            baflist = []
            bin = bins.readline().strip().split("\t")
            [chrom, start, end, na] = bin
        baflist.append(line[3])

while len(bin) > 1:
    if baflist:
        print("\t".join([chrom, start, end, ";".join(baflist)]))
        baflist = []
    else:
        print("\t".join([chrom, start, end, "NA"]))
    bin = bins.readline().strip().split("\t")
bins.close()
            
        
