#! /usr/bin/env Rscript
# prep_ploidetect2.R
#
# preprocesses input .bed file into an .RDS containing aggregated, pre-segmented data
#
# docopt script docstring
' prep_ploidetect2.R

Usage: 
prep_ploidetect2.R -i input -o output.txt [-c cyto] [-s single]

Options:
-i	--input		input	input .bed data file
-o	--output	output	output file
-c	--cyto	cytos	cytoband file
-s	--single	single  T/F whether it is singlesample or not
' -> doc
#
# load libraries
library(docopt)
library(devtools)
library(Ploidetect)
#
# Parse arguments
args = docopt(doc)
print(args)
#
# Read .bed table
all_data = read.table(args$input, header = F, sep = "\t", skip = 1, stringsAsFactors = F)
#
# cytoband
if("cyto" %in% names(args)){
    cytos = args$cyto
}else{
    cytos = F
}
#
# preprocess data
if("single" %in% names(args)){
    print(args$single)
    if(is.null(args$single)){
        single = F
    }else if (is.character(single)){
        single = as.logical(args$single)
    }
}else{
    single = F
}
print(single)
out = ploidetect_presegment(all_data, centromeres = cytos, singlesample = single)
#
# Save preprocessed data to output .RDS file
saveRDS(out, args$output)
