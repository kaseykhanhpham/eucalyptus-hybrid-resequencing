# R script to convert a BIM file to a genome annotation BED3 file (not the PLINK BED file)
# Usage: Rscript bim_to_bed.r [infile] [outfile]
# infile: the BIM file to be converted
# outfile: the name to export the resulting BED file to

infile_name <- commandArgs(trailingOnly = TRUE)[1]
outfile_name <- commandArgs(trailingOnly = TRUE)[2]
infile <- read.table(infile_name, header = FALSE, sep = "\t")
loc_ids <- infile$V2
chr_names <- unlist(strsplit(loc_ids, ":"))[1]
end_pos <- unlist(strsplit(loc_ids, ":"))[2]
end_pos <- as.numeric(end_pos)
start_pos <- end_pos - 1

bed_table <- cbind(chr_names, start_pos, end_pos)
write.table(bed_table, outfile_name, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")