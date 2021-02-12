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

library(docopt)
library(devtools)
#devtools::load_all("/projects/lculibrk_prj/CNV/Ploidetect-stable/Ploidetect-package/")
library(Ploidetect)

setDTthreads(1)

args = docopt(doc)

#source("R/ploidetect2.R")

in_rds = readRDS(args$input)
in_models = read.table(file = args$models, sep = " ", stringsAsFactors=F, header = T)

print(in_models)

print(in_rds$maxpeak)
print(in_models$tp)
print(in_models$ploidy)

if(ncol(in_models) == 1){
 in_models = data.frame("tp" = 1, "ploidy" = 2)
}

result = ploidetect_cna_sc(all_data = in_rds$all_data, segmented_data = in_rds$segmented_data, tp = in_models$tp[1], ploidy = in_models$ploidy[1], maxpeak=in_rds$maxpeak, verbose = T, max_iters = args$size)

pdf(args$plots)
result$cna_plots
dev.off()

write.table(file = args$output, x = result$cna_data, col.names = T, row.names=F)
write.table(file = paste0(gsub(".txt", "", args$output), "_condensed.txt"), x = result$segged_cna_data, col.names = T, row.names = F, quote = F, sep = "\t")
