# quick script to reformat output of vcftools FST windows to match egglib output
# Usage: Rscript reformat_fst.r [infile]

infile_name <- commandArgs(trailingOnly = TRUE)[1]
infile <- read.table(infile_name, header = TRUE)

colnames(infile) <- c("chr", "start", "end", "n_variants", "weighted_FST", "mean_FST")
