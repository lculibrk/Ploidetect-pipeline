' annotate.R

Usage: 
annotate.R -i cna_condensed.txt -a annotations.gtf -o output.bed

Options:
-i --input	cna	cna_condensed.txt
-a --annotations	annotations	annotations.gtf
-o --output output	output	output.bed
' -> doc

library(docopt)
library(data.table)

args = docopt(doc)
## Read gtf file, skip comment lines
annotations = fread(args$annotations, skip = 5)

## Assign names to gtf
names(annotations) = c("chr",
											 "annotation",
											 "type",
											 "pos",
											 "end",
											 "",
											 "strand",
											 ".",
											 "info")

## Define gtf parser
parse_gtf = function(gff, colnames){
	## Filter for gene annotations
	gff = gff[type == "gene"]
	## Get column names from info using row 1 as the example
	eg = gff$info[1]
	## Extract column names from leading word
	nms = gsub("([^ ]*) .*", "\\1", unlist(strsplit(eg, "; ")))
	## Convert string to multiple columns
	gff[,paste0(nms):=tstrsplit(info, "; ")]
	## Extract strings within brackets, applied to the new columns
	gff[, 9:ncol(gff)] = gff[,lapply(.SD, function(x)gsub(".*\"(.*)\".*", "\\1", x)), .SDcols = 9:ncol(gff)]
	## Return desired columns
	return(gff[,..colnames])
}

## Parse gtf file
annotations = parse_gtf(annotations, colnames = c("chr", "pos", "end", "gene_name"))

## Read cnv file
cnvs = fread(args$input)

## Set key magic for data.table to do foverlaps magic
setkey(cnvs, chr, pos, end)
setkey(annotations, chr, pos, end)

## foverlaps magic
annotated = foverlaps(cnvs, annotations)

## Cleanup column names
names(annotated) = gsub("i\\.", "", names(annotated))

## Write output
fwrite(annotated, args$output, sep = "\t")