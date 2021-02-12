' run_ploidetect2.R

Usage: 
run_ploidetect2.R -i input -p plots.pdf -r rds.rds -o output.txt

Options:
-i --input input      input .RDS data file
-p --plots plots      output plots file
-o --output output    output file
-r --rds rds	      output .RDS
' -> doc

#if("docopt" %in% rownames(installed.packages())){
 library(docopt)
#}else{install.packages("docopt", repos='http://cran.us.r-project.org')}
#if("devtools" %in% rownames(installed.packages())){
 library(devtools)
#}else{install.packages("devtools", repos='http://cran.us.r-project.org')}
#devtools::load_all("Ploidetect/")

#if("Ploidetect" %in% rownames(installed.packages())){
 library(Ploidetect)
#}else{devtools::install_github("lculibrk/Ploidetect")}


args = docopt(doc)


print(paste0("run_ploidetect2.R -i ", args$input, " -p ", args$plots, " -r ", args$rds, " -o ", args$output))


#source("R/ploidetect2.R")

in_list = readRDS(args$input)

result = ploidetect(in_list)

write.table(result$model_info, args$output)

pdf(args$plots)
print(result$plots)
dev.off()

saveRDS(list("maxpeak" = result$maxpeak, "all_data" = in_list$all_data, "segmented_data" = result$segmented_data), args$rds)

warnings()
