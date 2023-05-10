# quick script to reformat output of vcftools FST windows to match egglib output
# Usage: Rscript reformat_fst.r [infile]

infile_name <- commandArgs(trailingOnly = TRUE)[1]
infile <- read.table(infile_name, header = TRUE)

colnames(infile) <- c("chr", "start", "end", "n_variants", "weighted_FST", "mean_FST")
infile$start <- as.numeric(infile$start - 1) # 0 indexing to match egglib

options(digits = 22)
write.table(infile, infile_name, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")