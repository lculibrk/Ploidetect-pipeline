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


#### Re-create CNV plots to allow for rasterization
#### to-do: implement this in the package (but I'm not touching the code for now, to avoid any rabbit holes)
CN_calls <- split(result$cna_data, f = result$cna_data)

# Check if human genome is being used compared to e.g. mm10 and sort numerically instead of lexicographically
if("22" %in% names(CN_calls)){
	ordering_vec = c(1:22, "X")
	CN_calls = CN_calls[ordering_vec]
}

  CN_palette <- c("0" = "#cc0000", 
                  "1" = "#000066", 
                  "2" = "#26d953", 
                  "3" = "#609f70", 
                  "4" = "#cccc00", 
                  "5" = "#80804d",
                  "6" = "#cc6600", 
                  "7" = "#856647", 
                  "8" = "#cc0000"
  )
  
CNA_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)[1]] | pos > centromeres$end[which(centromeres$chr %in% chr)[2]]) %>% ggplot(aes(x = pos, y = log(corrected_depth, base = 2), color = as.character(state))) + 
      geom_point_rast(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("log(Read Depth)") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " copy number profile")) + 
      theme_bw()
  })
vaf_plot <- lapply(CN_calls, function(x){
    chr = x$chr[1]
    x %>% filter(end < centromeres$pos[which(centromeres$chr %in% chr)][1] | pos > centromeres$end[which(centromeres$chr %in% chr)][2]) %>% filter(!is.na(maf)) %>% ggplot(aes(x = pos, y = unlist(unmerge_mafs_grouped(maf, flip = T)), color = as.character(state))) + 
      geom_point_rast(size = 0.5) + 
      scale_color_manual(name = "State",
                         values = CN_palette, 
                         labels = c("0" = "HOMD", 
                                    "1" = "CN = 1", 
                                    "2" = "CN = 2 HET", 
                                    "3" = "CN = 2 HOM", 
                                    "4" = "CN = 3 HET", 
                                    "5" = "CN = 3 HOM", 
                                    "6" = "CN = 4 HET", 
                                    "7" = "CN = 4 HOM", 
                                    "8" = "CN = 5+")) + 
      ylab("Major allele frequency") + 
      xlab("position") + 
      ggtitle(paste0("Chromosome ", chr, " allele frequency profile")) + 
      scale_y_continuous(limits = c(0.5, 1)) +
      theme_bw()
})
  
cna_plots <- list()
  
for(i in 1:length(CNA_plot)){
  cna_plots[i] <- list(plot_grid(CNA_plot[[i]], vaf_plot[[i]], align = "v", axis = "l", ncol = 1))
}
  
CN_calls <- do.call(rbind.data.frame, CN_calls)

## End to-do area
#
# Open pdf device
pdf(args$plots)
# print plots to send them to the pdf device
CN_calls
# Close pdf device
dev.off()
# Write CNV per-bin calls to file
write.table(file = args$output, x = CN_calls, col.names = T, row.names=F, sep = "\t")
# Write segmented CNV calls to file
write.table(file = paste0(gsub(".txt", "", args$output), "_condensed.txt"), x = result$segged_cna_data, col.names = T, row.names = F, quote = F, sep = "\t")


