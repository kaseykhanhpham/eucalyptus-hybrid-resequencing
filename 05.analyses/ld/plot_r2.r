# R script to plot output of average plink r2 values by distance
# takes a CSV file with distance between SNPs and their average R2
# usage: Rscript plot_r2.r [infile] [outfile] [chromosome name] [dist interval]

# import packages
library(ggplot2)

# initialize location variables
infile <- commandArgs(trailingOnly = TRUE)[1]
outfile <- commandArgs(trailingOnly = TRUE)[2]
chr_name <- commandArgs(trailingOnly = TRUE)[3]

# read input
r2_tab <- read.csv(infile, header = TRUE)

# plot
r2_plot <- ggplot(dat = r2_tab, aes(x = dist, y = r2)) + geom_line() + ggtitle(label = paste("Linkage Disequilibrium of ", chr_name, sep = ""), subtitle = "R2 vs distance (bp)") + xlab("distance (bp)") + ylab("R2")

png(filename = outfile, width = 1600, height = 1000, units = "px")
r2_plot
dev.off()
