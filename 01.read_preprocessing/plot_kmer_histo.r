# R script to automatically plot k-mer count results
# Author: Kasey Pham
# Date: 10/5/21
# based on https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html

# Usage: Rscript plot_kmer_histo.r [infile_name] [outfile_name] [start_index]

infile_name <- commandArgs(trailingOnly=TRUE)[1]
outfile_name <- commandArgs(trailingOnly=TRUE)[2]
start_index <- commandArgs(trailingOnly=TRUE)[3]

kmer_table <- read.table(infile_name)

png(outfile_name)
plot(kmer_table[start_index:200,], type = "l", main = paste())
dev.off()
