# R script to convert posterior table from Ancestry HMM into a BED
# genome annotation file
# Usage: Rscript post_to_bed.r [infile] [outfile]

options(scipen = 999)

infile_name <- commandArgs(trailingOnly = TRUE)[1]
outfile_name <- commandArgs(trailingOnly = TRUE)[2]

infile <- read.table(infile_name, header = TRUE, as.is = TRUE)
end_pos <- as.numeric(infile$position)
start_pos <- end_pos - 1

bed_tab <- cbind(infile$chrom, start_pos, infile$position)

write.table(bed_tab, outfile_name, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
