# R script to plot output of average plink r2 values by distance
# takes a CSV file with distance between SNPs and their average R2
# usage: Rscript plot_r2.r [infile] [outfile] [dist interval]

# import packages
library(ggplot2)

# initialize location variables
infile <- commandArgs(trailingOnly = TRUE)[1]
outfile <- commandArgs(trailingOnly = TRUE)[2]
zoom <- FALSE

# check for optional x-axis interval specification and get interval
if(length(commandArgs(trailingOnly = TRUE)) == 3){
    dist_interval <- commandArgs(trailingOnly = TRUE)[3]
    start_int <- as.numeric(unlist(strsplit(dist_interval, ","))[1])
    end_int <- as.numeric(unlist(strsplit(dist_interval, ","))[2])
    zoom <- TRUE
}

# read input
r2_tab <- read.csv(infile, header = TRUE)

# plot
r2_plot <- ggplot(dat = r2_tab, aes(x = dist, y = r2)) + geom_line() + theme(text=element_text(size = 30))
r2_plot <- r2_plot + xlab("distance (bp)") + ylab("R2")

if(zoom){
    r2_plot <- r2_plot + ggtitle(label = paste("Linkage Disequilibrium", sep = ""), subtitle = paste("Zoomed - ", start_int, "bp to ", end_int, "bp", sep = ""))
    r2_plot <- r2_plot + xlim(start_int, end_int)
} else {
    r2_plot <- r2_plot + ggtitle(label = paste("Linkage Disequilibrium", sep = ""), subtitle = "Full view")
}

png(filename = outfile, width = 1600, height = 1000, units = "px")
r2_plot
dev.off()
