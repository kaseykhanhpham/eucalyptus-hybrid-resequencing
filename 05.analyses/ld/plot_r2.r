# R script to plot output of average plink r2 values by distance
# takes a CSV file with distance between SNPs and their average R2
# usage: Rscript plot_r2.r [infile] [outfile] [chromosome name] [dist interval]

# import packages
library(ggplot2)

# initialize location variables
infile_dir <- commandArgs(trailingOnly = TRUE)[1]
outfile_dir <- commandArgs(trailingOnly = TRUE)[2]
chr_name <- commandArgs(trailingOnly = TRUE)[3]
dist_int <- commandArgs(trailingOnly = TRUE)[4]

setwd(working_dir)

# read input
r2_tab <- read.csv(infile_dir, header = TRUE)

avg_r2 <- sapply(seq(from = 1, to = nrow(r2_tab), by = dist_int), function(x) mean(r2_tab[seq(from = x, to = (x + (dist_int - 1))), "r2"]))

avg_windows <- as.data.frame(cbind(dist = seq(from = dist_int, to = nrow(r2_tab), by = dist_int), avg_r2))

# plot
r2_plot <- ggplot(dat = avg_windows, aes(x = dist, y = avg_r2)) + geom_line() + ggtitle(label = paste("Linkage Disequilibrium of ", chr_name, sep = ""), subtitle = "R2 vs distance (bp)") + xlab("distance (bp)") + ylab("R2")

png(width = 1600, height = 1000)
r2_plot
dev.off()
