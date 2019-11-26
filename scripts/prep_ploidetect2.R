' prep_ploidetect2.R

Usage: 
prep_ploidetect2.R -i input -o output.txt

Options:
-i --input input      input .bed data file
-o --output output    output file
' -> doc
library(docopt)
library(devtools)
library(Ploidetect)
args = docopt(doc)

all_data = read.table()
out = ploidetect_presegment(all_data)
saveRDS(out, args$output)

#read.table(args$input, header = F, sep = "\t", skip = 1, stringsAsFactors = F)