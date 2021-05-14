#! /usr/bin/env Rscript
# ploidetect_copynumber.R
#
# Runs Ploidetect's copy number caller on Ploidetect-preprocessed input data using Ploidetect's tumor purity/ploidy estimates

# Docopt script docstring
' ploidetect_copynumber.R

Usage:
ploidetect_copynumber.R -i input -m models.txt -p plots.pdf -o output.txt [--size=SIZE]

Options:
-i	--input	input	input .RDS data file
-m	--models	models	input .txt file of desired models for copy number calling
-p	--plots	plots	output plots file
-o	--output	output	output file
		--size	maximum iterations to increase resolution [default: Inf]

' -> doc
#
# Load packages
library(docopt)
library(devtools)
library(data.table)
library(Ploidetect)
library(ggrastr)
#
# Force data.table to use only one thread (defaults to a large number)
setDTthreads(1)
#
# Parse arguments
args = docopt(doc)
#
# Read preprocessed input
in_rds = readRDS(args$input)
#
# Read purity/ploidy models
in_models = read.table(file = args$models, sep = " ", stringsAsFactors=F, header = T)
#
# Load centromere positions packaged with Ploidetect
data(centromeres)
#
# Check if Ploidetect couldn't detect purity/ploidy from cnv information and default to 100% purity and diploid
if(ncol(in_models) == 1){
 in_models = data.frame("tp" = 1, "ploidy" = 2)
}
#
# Run ploidetect_cna_sc (subclone aware CNV caller)
result = ploidetect_cna_sc(all_data = in_rds$all_data, segmented_data = in_rds$segmented_data, tp = in_models$tp[1], ploidy = in_models$ploidy[1], maxpeak=in_rds$maxpeak, verbose = T, max_iters = args$size)
print(str(result))
print(names(result))
# Get the cna plots from result object
cna_plots = result$cna_plots
# Get CNV objects from result object
CN_calls = result$cna_data
# Open pdf device
pdf(args$plots)
# print plots to send them to the pdf device
cna_plots
# Close pdf device
dev.off()
# Write CNV per-bin calls to file
write.table(file = args$output, x = CN_calls, col.names = T, row.names=F, sep = "\t")
# Write segmented CNV calls to file
write.table(file = paste0(gsub(".txt", "", args$output), "_condensed.txt"), x = result$segged_cna_data, col.names = T, row.names = F, quote = F, sep = "\t")
