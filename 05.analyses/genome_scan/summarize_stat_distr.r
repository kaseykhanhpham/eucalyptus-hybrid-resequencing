# R script to summarize the output of sliding window calculations of stats in the python package egglib
# Usage: Rscript summarize_stat_distr.r [infile] [stat] [minimum snps]
# infile: table output of egglib or vcftools to be processed for outlier windows
# stat: the stat of interest, should match the column reporting the stat for the window
# min_snps: The minimum number of SNPs to include a window in calculations

infile_name <- commandArgs(trailingOnly = TRUE)[1]
stat <- commandArgs(trailingOnly = TRUE)[2]
min_snps <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

# import input and filter by num SNPs
in_table <- read.table(infile_name, header = TRUE, sep = "\t", na.strings="")
filtered_table <- in_table[which(in_table$num_var >= min_snps),]

# Report general stats to stdout
write(paste("File:", infile_name, ":"), stdout())
write("min\t1stquart\tmedian\tmean\t3rdquart\tmax", stdout())
write(summary(filtered_table[, stat]), stdout())
write(paste("sd:", SD = sd(filtered_table[, stat])), stdout())
write("", stdout())
# Print stat distribution to PNG file
png(paste(infile_name, "_distr.png", sep = ""), width = 500, height = 500)
    hist(filtered_table[, stat], main = paste(stat, "distribution"), xlab = stat)
dev.off()