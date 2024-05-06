# Rscript to plot windows of interest across multiple genome scan metrics
# for introgression, selection, etc.
# Usage: Rscript plot_chr_wins.r -c [chr] -s [chr_size] -o [out_prefix] -y [dxy] -d [dsuite] -a [ahmm_list] -e [elai_dose] -i [elai_snpinfo] -g [elai_fam]
# Last updated: 3/26/2024

# Import libs
options(scipen = 999) # don't use scientific notation
library(argparser)
win_funs_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/scripts" # CHANGE THIS TO WHERE YOU ARE KEEPING THE FUNCTIONS FILE
win_funs <- paste(win_funs_dir, "win_processing_funs.r", sep = "/")
win_funs_2 <- paste(win_funs_dir, "intr_lais_funs.r", sep = "/")
source(win_funs)
source(win_funs_2)

# Parse arguments
parser <- arg_parser("A script to plot genome scan windows of interest for each chromosome.")

parser <- add_argument(parser, "-c", help = "Chromosome name", default = "Chr01", type = "character")
parser <- add_argument(parser, "-s", help = "Chromosome size", default = 42219553, type = "numeric")
parser <- add_argument(parser, "-o", help = "Output prefix", default = "Chr01", type = "character")
parser <- add_argument(parser, "-y", help = "Dxy filename", default = "all_dxy.txt", type = "character")
parser <- add_argument(parser, "-d", help = "Dsuite filename", default = "local_Fstats.txt", type = "character")
parser <- add_argument(parser, "-a", help = "List of AHMM files", default = "ahmm_list.txt", type = "character")
parser <- add_argument(parser, "-e", help = "ELAI dose file", default = "Chr01_avg_ps12.txt", type = "character")
parser <- add_argument(parser, "-i", help = "ELAI SNP info file", default = "Chr01_r1_snpinfo.txt", type = "character")
parser <- add_argument(parser, "-g", help = ".fam file used to generate ELAI input for introgressed samples", default = "samples.fam", type = "character")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))
chr <- args_list[["c"]]
chr_size <- args_list[["s"]]
out_pref <- args_list[["o"]]
dxy_name <- args_list[["y"]]
dsuite_name <- args_list[["d"]]
ahmm_list_name <- args_list[["a"]]
elai_dose_name <- args_list[["e"]]
elai_snpinfo_name <- args_list[["i"]]
elai_samples_name <- args_list[["g"]]

# import AHMM file list
ahmm_list <- read.table(ahmm_list_name, header = FALSE)$V1
# import ELAI samples
elai_samples <- read.table(elai_samples_name, header = FALSE)$V1

# Make chromosome plots
# NOTE: LEGENDS AREN'T CURRENTLY GENERATED, YOU'LL HAVE TO ADD THAT SEPARATELY AFTER COMBINING THE PLOTS.

# Ancestry_HMM introgression
message("Plotting dxy peaks and df...", stderr())
png(paste(out_pref, chr, "intr_alt.png", sep = "_"), width = 1000, height = 3500)
    plot_intr_alt(ahmm_tabnames = ahmm_list, ahmm_outname = paste(out_pref, chr, "ahmm.bed", sep = "_"), dxy_tabname = dxy_name, dsuite_tabname = dsuite_name, scan_outname = paste(out_pref, chr, "gscan_alt.bed", sep = "_"), chr = chr, chr_size = chr_size)
dev.off()
