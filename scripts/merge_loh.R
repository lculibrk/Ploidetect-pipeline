' merge_loh.R

Usage:
merge_cha.R -l loh -w windows -o output

Options:
-l --loh loh          input LOH file
-w --windows windows  input window file
-o --output output    output file
' -> doc

library(docopt)
library(dplyr)
library(GenomicRanges)

# Use script-relative import replacement for
#   source("scripts/merge_mafs.R")
this_folder <- funr::get_script_path()
source(paste(this_folder, "merge_mafs.R", sep="/"))

args <- docopt(doc)
print(args)
if(args$loh == "STDIN"){
  print("Reading from STDIN")
  loh <- read.table(file("stdin"), header=F, sep = "\t", stringsAsFactors = F)
}else{loh <- read.table(args$loh, header=F, sep="\t", stringsAsFactors = F)}

windows <- read.table(args$windows, header = F, sep = "\t", stringsAsFactors = F)

## Change loh and windows to GRanges format
loh <- GRanges(seqnames = loh$V1, ranges = IRanges(start = loh$V2+1, end = loh$V3), vaf = loh$V4)
windows <- GRanges(seqnames=windows$V1, ranges = IRanges(start = windows$V2+1, end = windows$V3))
## Find overlaps between the two datasets
overlaps <- findOverlaps(windows, loh)
## Note the rows
loh$rows <- 1:length(loh)
windows$rows <- 1:length(windows)
## Convert to dataframes
loh <- as.data.frame(loh)
windows <- as.data.frame(windows)
## Create a mapping df
mapping <- data.frame("mapping" = overlaps@from, "rows" = overlaps@to)
## Join mapping with loh
loh <- loh %>% left_join(mapping, by = "rows")
loh <- as.data.frame(loh) %>% dplyr::select(seqnames, start, end, vaf, mapping) %>% rename("seqnames" = "chr", "start" = "pos")
## Aggregate by window
loh <- loh %>% group_by(chr, mapping) %>% summarise(pos = min(pos) - 1, end = max(end), vaf = merge_mafs(vaf, na.rm = T, exp = T))
## Map back to windows
windows <- windows %>% left_join(loh[,c(2, 5)], by = c("rows" = "mapping")) %>% rename("seqnames" = "chr", "start" = "pos") %>% dplyr::select(chr, pos, end, vaf) %>% mutate(pos = pos - 1)
## Write to file
write.table(windows, args$output, col.names = F, row.names = F, quote = F, sep = "\t")
