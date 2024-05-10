#!/usr/bin/env Rscript

"
Purpose: R script to plot coverage from a table of coverage in genome-wide
sliding windows.
Input: A list of input tables of average coverage across windows for BAM files
       Input table should be tab delimited and have position in column 1 and 
       average coverage in column 2
Usage: plot_coverage.r -i [infile_list] -n [sample_table] -o [outfile]
                       -m [chrnames] -p [chrsizes] -s [winsize] -c [offset]
Dependencies: argparser, ggplot2
"

library(argparser)
library(ggplot2)

parser <- arg_parser("A script to plot sequencing coverage")

parser <- add_argument(parser, "-i", help = "list of avg coverage input files",
                       default = "infiles.txt", type = "character")
parser <- add_argument(parser, "-n", help = "table of sample names and pop'n",
                       default = "names_family.tab", type = "character")
parser <- add_argument(parser, "-o", help = "output name in PNG format",
                       default = "cover.png", type = "character")
parser <- add_argument(parser, "-m", help = "list of chromosome names",
                       default = "chr_names.txt", type = "character")
parser <- add_argument(parser, "-p", 
                       help = "list of chromosome sizes",
                       default = "chr_sizes.txt", type = "character")
parser <- add_argument(parser, "-s", help = "window size (bp)",
                       default = 100000, type = "numeric")
parser <- add_argument(parser, "-c", help = "offset size (bp)",
                       default = 100000000, type = "numeric")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

## IMPORT AND FORMAT SAMPLE TABLES
# import file lists to read
cov_tab_name_list <- read.table(args_list[["i"]], header = FALSE,
                     as.is = TRUE)$V1
sample_name_tab <- read.table(args_list[["n"]], header = TRUE,
                              as.is = TRUE)
"
FILE FORMAT NOTE:
sample population table should be tab delim with a header, eg:
name    pop
sample1 pop1
sample2 pop2
sample3 pop1
sample4 pop1
sample5 pop2
"

# import and sort input tables by base pair position
cov_tab_list <- list()
for (i in seq_len(length(cov_tab_name_list))) {
  samp <- sample_name_tab[i, "name"]
  # import coverage table for sample
  in_tab <- read.table(cov_tab_name_list[i], header = FALSE,
                       col.names = c("pos", "cov"))
  # sort table by window position (first column)
  cov_tab_list[[samp]] <- in_tab[order(in_tab$pos), ]
}

## CREATE WINDOW BINS FOR AVERAGING
# import input information
chr_list <- read.table(args_list[["m"]], header = FALSE, as.is = TRUE)$V1
chr_sizes <- as.numeric(read.table(args_list[["p"]], header = FALSE)$V1)
winsize <- args_list[["s"]]
offs <- args_list[["c"]]

# initialize storage vars
window_starts <- c()
window_ends <- c()
chr_centers <- c()

# loop through chr to populate storage vars with window starts and ends
for (i in seq_len(length(chr_sizes))) {
  # calculate starts by adding num windows tiling chr size to offset
  num_chr_starts_aftr <- floor(chr_sizes[i] / winsize)
  first_s <- offs * (i - 1)
  chr_starts <- c(first_s, (first_s + (winsize * seq_len(num_chr_starts_aftr))))
  # calculate ends as 1 less than start of next window
  chr_ends <- chr_starts[2:length(chr_starts)] - 1
  chr_ends <- c(chr_ends, chr_sizes[i])
  # append to storage vars
  window_starts <- c(window_starts, chr_starts)
  window_ends <- c(window_ends, chr_ends)
  # calculate center as half of chr size + offset
  chr_centers <- c(chr_centers, round(chr_starts[1] + (chr_sizes[i] / 2)))
}

## AVERAGE COVERAGE FOR A CERTAIN WINDOW ACROSS SAMPLES IN A POP'N
pop_list <- unique(sample_name_tab)
# initiate storage vars for each pop avg table
all_pos <- c()
all_cov <- c()
all_pop <- c()

# do averaging
for (pop in pop_list) {
  samps <- sample_name_tab[which(sample_name_tab$pop == pop), "name"]
  # initialize storage var for coverages
  pos <- c()
  cov <- c()
  # across each window in genome, average over all samples in the pop
  for (i in seq_len(length(window_starts))) {
    curr_win_vals <- c() # storage for window avg across samples
    for (s in samps) {
      # entries bigger than win start
      starts_after <- which(cov_tab_list[[s]]$pos >= window_starts[i])
      # entries smaller than win end
      ends_before <- which(cov_tab_list[[s]]$pos <= window_ends[i])
      # entries flagged for both to be considered for average
      in_win <- intersect(starts_after, ends_before)
      # add valid coverages to storage var if there are any for sample
      if (length(in_win) > 0) {
        curr_win_vals <- c(curr_win_vals, cov_tab_list[[s]][in_win, "cov"])
      }
    }
    # add position and calculate mean if any samples had entries in the win
    if (length(curr_win_vals) > 0) {
      pos <- c(pos, window_starts[i])
      cov <- c(cov, mean(curr_win_vals))
    }
  }
  # add population means for wins to final avg storage vars
  all_pos <- c(all_pos, pos)
  all_cov <- c(all_cov, cov)
  all_pop <- c(all_pop, rep(pop, length(pos)))
}
# build pop avg table (for plotting in ggplot2)
pop_avg_tab <- data.frame(pos = all_pos, cov = all_cov, pop = all_pop)

## PLOT COVERAGE PER POPULATION PER WINDOW
plot_title <- "Average Coverage per Window per Population"

png(args_list[["o"]], height = 1600, width = 1800, units = "px")
# scatterplot of true coverage values
# NOTE: scatterplot is restricted to 200x coverage on the y-axis.
# this will exclude sites with higher coverage from visualization!

# define colors to use for each population
colscale <- c("goldenrod1", "#606060", "deepskyblue4")
# plot chromosome names and positions if provided
plt <- ggplot(pop_avg_tab, aes(x = pos, y = cov, col = pop))
plt <- plt + theme_bw(base_size = 40)
plt <- plt + geom_point(alpha = 0.8, size = 2) # transparency 80%
plt <- plt + scale_colour_manual(values = colscale, name = "Population",
                                 labels = c("E. cordata", "MR E. globulus",
                                            "Pure E. globulus"))
plt <- plt + ggtitle(plot_title)
plt <- plt + labs(x = "Window Position", y = "Average Coverage per Window")
plt <- plt + ylim(0, 200) # RESTRICT TO PLOTTING UP TO 200X COVERAGE
# set x-axis breakpoints to chromosomes
plt <- plt + scale_x_continuous(breaks = chr_centers, labels = chr_list)
# add line for 30x coverage
plt <- plt + geom_hline(yintercept = 30, linetype = "dashed",
                        color = "#000000")

print(plt)
dev.off()