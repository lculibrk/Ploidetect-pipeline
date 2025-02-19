#! /usr/bin/env Rscript
# prep_ploidetect2.R
#
# preprocesses input .bed file into an .RDS containing aggregated, pre-segmented data
#
# docopt script docstring
' prep_ploidetect2.R

Usage: 
prep_ploidetect2.R -i input -o output.txt [-c cyto]

Options:
-i	--input		input	input .bed data file
-o	--output	output	output file
-c	--cyto	cytos	cytoband file
' -> doc
#
# load libraries
library(docopt)
library(devtools)
library(Ploidetect)
#
# Parse arguments
args = docopt(doc)

print(paste0(c(args$input,
              " sent to Ploidetect::ploidetect_presegment version: ",
              packageVersion("Ploidetect"))))
# Read .bed table
all_data = read.table(args$input, header = F, sep = "\t", skip = 1, stringsAsFactors = F)
#
# cytoband
if("cyto" %in% names(args)){
    cytos = args$cyto
}else{
    cytos = F
}
print(cytos)
#
# preprocess data
out = ploidetect_presegment(all_data, centromeres = cytos)
#
# Save preprocessed data to output .RDS file
saveRDS(out, args$output)
