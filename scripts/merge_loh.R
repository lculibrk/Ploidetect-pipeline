# merge_loh.R
#
# Inputs:  a .bed-formatted file of BAFs by position (only stdin has been 
#          tested) and a windows file describing the bins to map to
# Outputs: a .bed-formatted file with aggregated BAFs for each provided bin
# 
# Performs BAF flipping and aggregation. 
# Since BAFs are symmetrical about 0.5, attempting to take the mean or even
# median can erroneously return 0.5. This script counts the sub-0.5 and 
# super-0.5 BAF values, and flips the category with fewer BAFs about 0.5
# and takes the mean value
#
# TODO: Replace dplyr with data.table
#
#
# Docopt doc for interpreting script invocation
' merge_loh.R

Usage:
merge_cha.R -l loh -w windows -o output

Options:
-l --loh loh          input LOH file
-w --windows windows  input window file
-o --output output    output file
' -> doc
#
# Load libraries
library(docopt)
library(dplyr)
library(GenomicRanges)
#
# Scipen to prevent scientific notation
options(scipen=999)
#
# Read arguments
args <- docopt(doc)
#
# Read input either from stdin or file
# Check if input is stdin
if(args$loh == "STDIN"){
  # Debugging message
  print("Reading from STDIN")
  # Read file from stdin
  loh <- read.table(file("stdin"), header=F, sep = "\t", stringsAsFactors = F)
  # Otherwise, read from provided file path
}else{loh <- read.table(args$loh, header=F, sep="\t", stringsAsFactors = F)}
#
# merge_mafs.R functions
merge_mafs <- function(inputmafs, na.rm = T, flip = F, exp = F){
  # Takes in vector of input BAF valus and flips them about 0.5 according to
  # majority vote between sub-0.5 and super-0.5
  #
  # Deprecated functionality, should test & remove at some point
  if(exp){
    if(length(na.omit(inputmafs)) == 0){
      return(NA_character_)
    }
    return(paste(na.omit(inputmafs), collapse = ";"))
  }
  # Remove NA values
  if(na.rm){
    inputmafs <- inputmafs[which(!is.na(inputmafs))]
  }
  # Center values at 0.5
  mafs_about_zero <- inputmafs - 0.5
  # Count how many are sub-0.5 and super-0.5
  negatives <- which(mafs_about_zero < 0)
  positives <- which(mafs_about_zero > 0)
  # Compute absolute magnitude about 0.5
  flipped_mafs <- abs(mafs_about_zero)
  # If most of the BAFs are >0.5, flip to >0.5 and obtain mean
  if(length(positives) > length(negatives)){
    out <- mean(0.5 + flipped_mafs)
  # Otherwise flip to <0.5 and obtain mean
  }else{
    out <- mean(0.5 - flipped_mafs)
  }
  # Flip if specified
  if(flip){
    out <- abs(out - 0.5) + 0.5
  }
  return(out)
}
#
# Unused
unmerge_mafs <- function(merged_mafs, flip = F){
  mafs <- na.omit(merged_mafs)
  if(flip){
    return(abs(as.numeric(unlist(lapply(mafs, strsplit, split = ";"), recursive = T)) - 0.5) + 0.5)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}
#
# Unused
flip_merged_mafs <- function(merged_mafs){
  mafs <- unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
  flipped_mafs <- lapply(mafs, function(x)merge_mafs(abs(as.numeric(x) - 0.5) + 0.5, exp = T))
  return(flipped_mafs)
}
#
# Unused
unmerge_mafs_grouped <- function(merged_mafs, flip = F){
  if(flip){
    mafs <- unlist(lapply(merged_mafs, strsplit, split = ";"), recursive = F)
    mafs <- lapply(mafs, function(x){median(abs(as.numeric(x) - 0.5) + 0.5)})
    return(mafs)
  }
  unlist(lapply(mafs, strsplit, split = ";"), recursive = F)
}
#
# Read bins to map to
windows <- read.table(args$windows, header = F, sep = "\t", stringsAsFactors = F)
#
# Change loh and bins to GRanges format
loh <- GRanges(seqnames = loh$V1, ranges = IRanges(start = loh$V2+1, end = loh$V3), vaf = loh$V4)
windows <- GRanges(seqnames=windows$V1, ranges = IRanges(start = windows$V2+1, end = windows$V3))
# Find overlaps between the two datasets
overlaps <- findOverlaps(windows, loh)
# Note the rows
loh$rows <- 1:length(loh)
windows$rows <- 1:length(windows)
# Convert to dataframes
loh <- as.data.frame(loh)
windows <- as.data.frame(windows)
# Create a mapping df
mapping <- data.frame("mapping" = overlaps@from, "rows" = overlaps@to)
# Join mapping with loh
loh <- loh %>% left_join(mapping, by = "rows")
loh <- as.data.frame(loh) %>% dplyr::select(seqnames, start, end, vaf, mapping) %>% rename("seqnames" = "chr", "start" = "pos")
# Aggregate by window
loh <- loh %>% group_by(chr, mapping) %>% summarise(pos = min(pos) - 1, end = max(end), vaf = merge_mafs(vaf, na.rm = T, exp = T))
# Map back to windows
windows <- windows %>% left_join(loh[,c(2, 5)], by = c("rows" = "mapping")) %>% rename("seqnames" = "chr", "start" = "pos") %>% dplyr::select(chr, pos, end, vaf) %>% mutate(pos = pos - 1)
# Write to file
write.table(windows, args$output, col.names = F, row.names = F, quote = F, sep = "\t")
