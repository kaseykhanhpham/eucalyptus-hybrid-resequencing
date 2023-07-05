# R script to retrieve window output of Ancestry_HMM heterozygous
# or homozygous for _E. cordata_ ancestry (1,1 or 0,2)
# Usage: Rscript get_cordata_posteriors.r [infile_name] [outfile_root] [cutoff]
#     infile_name: the table output of Ancestry_HMM for an admixed individual
#     outfile_root: the base for the output files' names
#     cutoff: a double, the minimum cutoff posterior probability for a site to be included

infile_name <- commandArgs(trailingOnly = TRUE)[1]
outfile_root <- commandArgs(trailingOnly = TRUE)[2]
cutoff <- as.numeric(commandArgs(trailingOnly = TRUE)[3])

infile <- read.table(infile_name, header = TRUE, sep = "\t")
heterozygotes <- infile[which(infile[,4] > cutoff),]
homozygotes <- infile[which(infile[,5] > cutoff),]

het_outname <- paste(outfile_root, "_het_", cutoff, ".tab", sep = "")
hom_outname <- paste(outfile_root, "_hom_", cutoff, ".tab", sep = "")

write.table(heterozygotes, het_outname, quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(homozygotes, hom_outname, quote = FALSE, row.names = FALSE, col.names = TRUE)
