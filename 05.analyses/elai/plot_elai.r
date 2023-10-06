# R script to plot the results of ELAI for an individual
# Usage: R plot_elai.r [output_dir] [prefix_pattern] [chr_list] [ind_list] [ind] [output_name]
#     output_dir: the output directory created by ELAI for the run
#     prefix_pattern: pattern for output file prefixes for the run with {CHR} standing in for the chromosome ID
#     chr_list: a file with one chromosome ID per line for the genome
#     ind_list: a list giving the order that individuals are in from left to right in the ELAI genotype input file.
#     ind: the individual ID to focus on for the plotting, should correspond to the ID in the genotype file.
#     output_name: name of the PDF file of graphs for each chromosome

library(ggplot2)
options(scipen = 999)

output_dir <- commandArgs(trailingOnly = TRUE)[1]
prefix_pattern <- commandArgs(trailingOnly = TRUE)[2]
chr_listname <- commandArgs(trailingOnly = TRUE)[3]
ind_listname <- commandArgs(trailingOnly = TRUE)[4]
ind <- commandArgs(trailingOnly = TRUE)[5]
outname <- commandArgs(trailingOnly = TRUE)[6]

# import input files from cmd line
chr_list <- read.table(chr_listname, header = FALSE, as.is = TRUE)$V1
ind_list <- read.table(ind_listname, header = FALSE, as.is = TRUE)$V1
indv_index <- which(ind_list == ind)

# open file for printing
write(paste("Plotting to", outname), stdout())
pdf(outname, onefile = TRUE)
# plot each chromosome separately
for(chr in chr_list){
    # import ELAI output files for the chromosome
    dosage_filename <- paste(output_dir, "/", gsub("{CHR}", chr, prefix_pattern, fixed = TRUE), ".ps21.txt", sep = "")
    dosage_in <- scan(dosage_filename)
    pos_filename <- paste(output_dir, "/", gsub("{CHR}", chr, prefix_pattern, fixed = TRUE), ".snpinfo.txt", sep = "")
    pos_in <- read.table(pos_filename, header = TRUE)

    dim(dosage_in) <- c(2,nrow(pos_in),20) # 2 ancestral populations, num snps, 20 inds

    # subset to just E. cordata allele dosages for the focal individual
    indv_doses <- dosage_in[2,,indv_index]

    dose_by_pos <- data.frame(pos = pos_in$pos, doses = indv_doses)
    cplot <- ggplot(dose_by_pos, aes(x = pos, y = doses)) + geom_line() + theme_light() + ggtitle(paste(chr, "ancestry dosages")) + xlab("Position (bp)") + ylab("Dosage (E. cordata alleles)") + ylim(0,2)
    print(cplot)
    write(paste("\tFinished plotting", chr), stdout())
}

dev.off()
write(paste("Done writing to", outname), stdout())
