#!/usr/bin/env Rscript
# R script to plot coverage from a table of coverage in genome-wide sliding windows.
# Input table should be tab delimited and have position in column 1 and average coverage in column 2
# Usage: plot_coverage.r -i [infile] -o [outfile] -n [samplename] -c [chrlist] -x [chrpos]
# Dependencies: argparser

library(argparser)
parser <- arg_parser("A script to plot sequencing coverage")

parser <- add_argument(parser, "-i", help = "tab-delim input table, with positions in column 1 and coverage in column 2, no header", default = "cover.txt", type = "character")
parser <- add_argument(parser, "-o", help = "output plot in PNG format", default = "cover.png", type = "character")
parser <- add_argument(parser, "-n", help = "name of sample", default = "Sample", type = "character")
parser <- add_argument(parser, "-c", help = "list of chromosome names", default = "No file name provided :V", type = "character")
parser <- add_argument(parser, "-p", help = "list of chromosome center positions in input table", default = "No file name provided :V", type = "character")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

# import input table
cov_tab <- read.table(args_list[["i"]], header = FALSE, col.names = c("pos", "cov"))
cov_tab <- cov_tab[order(cov_tab$pos),] # sort table by window position (first column)

# import chromosome lists if present
if(args_list[["c"]] != "No file name provided :V"){
    chr_list <- as.character(read.table(args_list[["c"]], header = FALSE)$V1)
    chr_pos <- as.numeric(read.table(args_list[["p"]], header = FALSE)$V1)
}

# plot coverage and smoothing line
plot_title <- paste(args_list[["n"]], "Coverage")
png(args_list[["o"]], height = 1200, width = 1200, units = "px")
# scatterplot of true coverage values
# IMPORTANT: scatterplot is restricted to 200x coverage on the y-axis -- this will exclude sites with higher coverage from visualization!
if(args_list[["c"]] != "No file name provided :V"){ # plot chromosome names and positions if provided
    plot(cov_tab$pos, cov_tab$cov, pch = 19, main = plot_title, cex.main = 5, xlab = "", ylab = "", cex.axis = 3, cex.lab = 2.5, col = "black", ylim = c(0,200), xaxt = "n")
    axis(side = 1, at = chr_pos, labels = chr_list, cex.axis = 2, las = 2)
} else { # plot raw position values if chr names and positions not provided
    plot(cov_tab$pos, cov_tab$cov, pch = 19, main = plot_title, cex.main = 5, xlab = "Position", ylab = "Coverage", cex.axis = 3, cex.lab = 2.5, col = "black", ylim = c(0,200))
}
# add coverage cutoff line at 30x
abline(h = 30, col = "#1790a5", lwd = 4)
dev.off()
