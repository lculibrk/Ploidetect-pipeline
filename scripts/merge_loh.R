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

args <- docopt(doc)
print(args)
if(args$loh == "STDIN"){
  print("Reading from STDIN")
  loh <- read.table(file("stdin"), header=F, sep = "\t", stringsAsFactors = F)
}else{loh <- read.table(args$loh, header=F, sep="\t", stringsAsFactors = F)}

# merge_mafs.R functions
merge_mafs <- function(inputmafs, na.rm = T, flip = F, exp = F){
  if(exp){
    if(length(na.omit(inputmafs)) == 0){
      return(NA_character_)
    }
    return(paste(na.omit(inputmafs), collapse = ";"))
  }
  if(na.rm){
    inputmafs <- inputmafs[which(!is.na(inputmafs))]
  }
  ## Compute magnitude about 0.5
  mafs_about_zero <- inputmafs - 0.5
  negatives <- which(mafs_about_zero < 0)
  positives <- which(mafs_about_zero > 0)
  ## Compute absolute magnitude about 0.5 and add 0.5
  flipped_mafs <- abs(mafs_about_zero)
  if(length(positives) > length(negatives)){
    out <- mean(0.5 + flipped_mafs)
  }else{
    out <- mean(0.5 - flipped_mafs)
  }
  if(flip){
    out <- abs(out - 0.5) + 0.5
  }
  return(out)
}
unmerge_mafs <- function(merged_mafs, flip = F){
  mafs <- na.omit(merged_mafs)
  if(flip){
    return(abs(as.numeric(unlist(lapply(mafs, strsplit, split = ";"), recursive = T)) - 0.5) + 0.5)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}
flip_merged_mafs <- function(merged_mafs){
  mafs <- unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
  flipped_mafs <- lapply(mafs, function(x)merge_mafs(abs(as.numeric(x) - 0.5) + 0.5, exp = T))
  return(flipped_mafs)
}

unmerge_mafs_grouped <- function(merged_mafs, flip = F){
  if(flip){
    mafs <- unlist(lapply(merged_mafs, strsplit, split = ";"), recursive = F)
    mafs <- lapply(mafs, function(x){median(abs(as.numeric(x) - 0.5) + 0.5)})
    return(mafs)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}

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
