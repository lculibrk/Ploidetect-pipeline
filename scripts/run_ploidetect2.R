# run_ploidetect2.R
#
# runs Ploidetect's tumor purity and ploidy caller on preprocessed output
#
# Takes input .RDS from prep_ploidetect2.R and performs tumor purity and ploidy
# calling. Outputs a tab-separated file containing purity/ploidy models and a 
# pdf containing a plot of the data distribution and a plot for each model
#
# docopt script docstring
' run_ploidetect2.R

Usage: 
run_ploidetect2.R -i input -p plots.pdf -r rds.rds -o output.txt [-c cyto]

Options:
-i --input input      input .RDS data file
-p --plots plots      output plots file
-o --output output    output file
-r --rds rds	      output .RDS
-c --cytos cytos      cytobands file
' -> doc
# load libraries
library(docopt)
library(devtools)
library(Ploidetect)
#
# Parse arguments
args = docopt(doc)
#
# show the script call as run
print(paste0("run_ploidetect2.R -i ", args$input, " -p ", args$plots, " -r ", args$rds, " -o ", args$output))
#
# Read input data .RDS file
in_list = readRDS(args$input)
#
# Run Ploidetect on the data
result = ploidetect(in_list)
#
# Write the output models to file
write.table(result$model_info, args$output, sep = "\t")
#
# Open pdf device
pdf(args$plots)
#
# Print plots to pdf device
print(result$plots)
#
# Close pdf device
dev.off()
#
# Write the data to .rds for debugging purposes
saveRDS(list("maxpeak" = result$maxpeak, "all_data" = in_list$all_data, "segmented_data" = result$segmented_data), args$rds)
#
# Print warnings
warnings()
