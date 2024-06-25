# Purpose: R script to automatically plot k-mer count results for each sample on
#          one plot.
# Author: Kasey Pham
# Revised: 2024/05/10
# based on https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html

# Usage: Rscript plot_kmer_histo.r -i [infile table] -o [outfile_name]
#        -s [start_index]

"Expected Infile Table Format:
file1   file2   name    pop
S1  S332    WA01    glob_MR
S2  S449    WA02    cord_MR
"

library(argparser)

parser <- arg_parser("A script to plot sequencing coverage")

parser <- add_argument(parser, "-i",
                       help = "table of input files and pop'n for each sample",
                       default = "sample_reads_table.txt", type = "character")
parser <- add_argument(parser, "-o", help = "name of output graph",
                       default = "kmers.png", type = "character")
parser <- add_argument(parser, "-s", help = "Starting x-axis value for plot",
                       default = 3, type = "numeric")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

# Import input table
intab <- read.table(args_list[["i"]], header = TRUE, as.is = TRUE, sep = "\t")

# Import kmer count files
infiles_list <- list()
for (i in seq_len(nrow(intab))) {
infiles_list[[intab[i, "name"]]][[1]] <- read.table(intab[i, "file1"])
infiles_list[[intab[i, "name"]]][[2]] <- read.table(intab[i, "file1"])
}

# define function for merging frequency counts from two files
# for one sample
sum_freqs <- function(list_item, pos){
  row_1 <- which(list_item[[1]]$V1 == pos)
  row_2 <- which(list_item[[2]]$V1 == pos)
  freq_sum <- sum(list_item[[1]]$V2[row_1],
                  list_item[[2]]$V2[row_2], na.rm = TRUE)
  return(freq_sum)
}

# populate final table for plotting each sample
kmer_tables_combined <- list()
max_freq <- 0
for (samp in names(infiles_list)) {
  # only consider k-mer sizes in both files for a single sample
  pos_combined <- order(intersect(infiles_list[[samp]][[1]]$V1,
                                  infiles_list[[samp]][[2]]$V2))
  # combine k-mer counts from both files for a single sample
  freq_combined <- unlist(sapply(pos_combined, function(p) 
                                               sum_freqs(infiles_list[[samp]],
                                               p)))
  # log-transform frequencies to plot diff scales together
  freq_combined_log <- log10(freq_combined)
  # record maximum frequency count for plotting
  if (max(freq_combined_log) > max_freq) {
    max_freq <- max(freq_combined_log)
  }
  kmer_tables_combined[[samp]] <- data.frame(pos = pos_combined,
                                             freq = freq_combined_log)
}

# log-transform frequencies for plotting


# plot k-mer frequency counts for each sample
tiff(args_list[["o"]], height = 6, width = 6, units = "in", res = 350)
# establish plot canvas for all samples with padding
plot(x = c(args_list[["s"]], 1000),
     y = c(0, round(max_freq + max_freq * 0.025)), col = "white",
     main = "Log10 Read k-mer Frequencies of All Samples",
     xlab = "K-mer Size", ylab = "log10(Frequency)",
     cex.main = 1.25, cex.axis = 1.25)
# loop through all samples and plot line with 80% opacity and line color by pop
for (samp in names(kmer_tables_combined)) {
  samp_pop <- intab[which(intab$name == samp), "pop"]
  # assign color based on population affiliation
  if (samp_pop == "glob_MR") {
    samp_col <- "#13bdd7"
  } else if (samp_pop == "glob_pure") {
    samp_col <- "#044075"
  } else {
    samp_col <- "#ffcb3d"
  }
  # add sample's k-mer distribution to plot
  lines(x = kmer_tables_combined[[samp]]$pos,
        y = kmer_tables_combined[[samp]]$freq,
        type = "l", col = samp_col, alpha = 0.75, lwd = 0.75)
}

dev.off()
