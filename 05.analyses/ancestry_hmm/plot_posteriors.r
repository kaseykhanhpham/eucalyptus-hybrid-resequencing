# R script to plot posteriors from Ancestry_HMM output for each chromosome
# Usage: Rscript plot_posteriors.r [posteriors_file] [colors (optional)]
#        posteriors_file: output file for an admixed individual from Ancestry_HMM
#        colors (optional): comma-separated list of colors to use for each ancestry state; must match the number of total states

library(ggplot2)
library(tidyr)

infile_name <- commandArgs(trailingOnly = TRUE)[1]
infile <- read.table(infile_name, header = TRUE, sep = "\t")
num_states <- ncol(infile) - 2
state_cols <- colnames(infile)[3:(2 + num_states)]

# Get colors
if(length(commandArgs(trailingOnly = TRUE)) > 1){
    colors_list <- unlist(strsplit(commandArgs(trailingOnly = TRUE)[2], ","))
} else {
    colors_list <- c()
}

chr_list <- unique(infile$chrom)

split_file_addr <- unlist(strsplit(infile_name, "/"))
infile_name_relative <- split_file_addr[length(split_file_addr)]
base_name <- gsub(".posterior$", "", infile_name_relative)
out_name <- paste(base_name, "posteriors.pdf", sep = "_")
write(paste("Plotting", infile_name_relative, "to", out_name), stdout())
pdf(file = out_name, onefile = TRUE)

# plot each chromosome separately
for(chr in chr_list){
    infile_subset <- infile[which(infile$chrom == chr), ]
    # merge state columns into one and add col for state
    subset_piv <- pivot_longer(infile_subset, cols = state_cols, names_to = "state", values_to = "posterior")
    # construct lineplot object in ggplot2
    cplot <- ggplot(subset_piv, aes(x = position, y = posterior)) + geom_line(aes(col = state)) + theme_light() + ggtitle(paste(chr, "ancestry posteriors"))
    # generate color palettes
    if(length(colors_list) > 0){
        cplot <- cplot + scale_color_manual(name = "State", labels = state_cols, values = colors_list)
    } else {
        cplot <- cplot + scale_colour_hue(name = "State", labels = state_cols)
    }
    # plot
    print(cplot)
    write(paste("\tFinished plotting", chr), stdout())
}
dev.off()

write(paste("Done writing to", out_name), stdout())
