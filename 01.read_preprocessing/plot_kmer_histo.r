# R script to automatically plot k-mer count results
# Author: Kasey Pham
# Date: 10/5/21
# based on https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html

# Usage: Rscript plot_kmer_histo.r [infile_name] [infile_name_2] [outfile_name] [start_index] [sample_name]

infile_1_name <- commandArgs(trailingOnly=TRUE)[1]
infile_2_name <- commandArgs(trailingOnly=TRUE)[2]
outfile_name <- commandArgs(trailingOnly=TRUE)[3]
start_index <- commandArgs(trailingOnly=TRUE)[4]
sample_name <- commandArgs(trailingOnly=TRUE)[5]

kmer_table_1 <- read.table(infile_1_name)
kmer_table_2 <- read.table(infile_2_name)
plot_title <- paste(sample_name, "read k-mer frequency")

pos_combined <- order(intersect(kmer_table_1$V1, kmer_table_2$V1))

sum_freqs <- function(pos){
    row_1 <- which(kmer_table_1$V1 == pos)
    row_2 <- which(kmer_table_2$V1 == pos)
    freq_sum <- sum(kmer_table_1$V2[row_1], kmer_table_2$V2[row_2], na.rm = TRUE)
    return(freq_sum)
}

freq_combined <- unlist(sapply(pos_combined, function(x) sum_freqs(x)))
kmer_table_combined <- data.frame(pos = pos_combined, freq = freq_combined)

png(outfile_name)
plot(kmer_table_combined[start_index:200,], type = "l", lwd = 2.5, main = plot_title, xlab = "", ylab = "", cex.main = 2.5, cex.axis = 2.5)
dev.off()
