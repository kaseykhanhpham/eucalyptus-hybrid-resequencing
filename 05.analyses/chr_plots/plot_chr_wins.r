"
Rscript to plot windows of interest across multiple genome scan metrics
for introgression, selection, etc.

Usage: Rscript plot_chr_wins.r -c [chr] -s [chr_size] -o [out_prefix]
-n [spp_diffs] -p [pi] -y [dxy] -f [fst] -l [ld] -r [recomb]
-t [tajimasd] -a [ahmm_list] -e [elai_dose] -i [elai_snpinfo] -g [elai_fam]

Last updated: 5/06/2024
"

# Import libs
options(scipen = 999) # don't use scientific notation
library(argparser)
# CHANGE VAR BELOW TO WHERE YOU ARE KEEPING THE FUNCTIONS FILE
win_funs_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
win_funs <- paste(win_funs_dir, "win_processing_funs.r", sep = "/")
source(win_funs)

# Define parsed arguments
parser <- arg_parser("A script to plot genome scan windows of interest for each
chromosome.")

parser <- add_argument(parser, "-c", help = "Chromosome name",
                       default = "Chr01", type = "character")
parser <- add_argument(parser, "-s", help = "Chromosome size",
                       default = 42219553, type = "numeric")
parser <- add_argument(parser, "-o", help = "Output prefix", default = "Chr01",
                       type = "character")
parser <- add_argument(parser, "-n", help = "Species differences filename",
                       default = "spp_diffs.tab", type = "character")
parser <- add_argument(parser, "-p", help = "Pi filename",
                       default = "all_pi.txt", type = "character")
parser <- add_argument(parser, "-y", help = "Dxy filename",
                       default = "all_dxy.txt", type = "character")
parser <- add_argument(parser, "-f", help = "Fst filename",
                       default = "all_fst.txt", type = "character")
# parser <- add_argument(parser, "-d", help = "Dsuite filename",
#                        default = "local_Fstats.txt", type = "character")
parser <- add_argument(parser, "-l", help = "LD filename",
                       default = "ld_windows.txt", type = "character")
parser <- add_argument(parser, "-r", help = "Recombination filename",
                       default = "recomb_windows.txt", type = "character")
parser <- add_argument(parser, "-t", help = "Tajima's D filename",
                       default = "snps.TajimaD", type = "character")
parser <- add_argument(parser, "-a", help = "List of AHMM files",
                       default = "ahmm_list.txt", type = "character")
parser <- add_argument(parser, "-e", help = "ELAI dose file",
                       default = "Chr01_avg_ps12.txt", type = "character")
parser <- add_argument(parser, "-i", help = "ELAI SNP info file",
                       default = "Chr01_r1_snpinfo.txt", type = "character")
parser <- add_argument(parser, "-g", help =
                       ".fam file used to make ELAI input for introgr samples",
                       default = "samples.fam", type = "character")

# retrieve arguments from parser
args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))
chr <- args_list[["c"]]
chr_size <- args_list[["s"]]
out_pref <- args_list[["o"]]
spp_diffs_name <- args_list[["n"]]
pi_name <- args_list[["p"]]
dxy_name <- args_list[["y"]]
fst_name <- args_list[["f"]]
# dsuite_name <- args_list[["d"]]
ld_name <- args_list[["l"]]
recomb_name <- args_list[["r"]]
tajd_name <- args_list[["t"]]
ahmm_list_name <- args_list[["a"]]
elai_dose_name <- args_list[["e"]]
elai_snpinfo_name <- args_list[["i"]]
elai_samples_name <- args_list[["g"]]

# import AHMM file list
ahmm_list <- read.table(ahmm_list_name, header = FALSE)$V1
# import ELAI samples [WAITING FOR ANALYSIS TO FINISH]
elai_samples <- read.table(elai_samples_name, header = FALSE)$V1

# Make chromosome plots
# NOTE: LEGENDS AREN'T CURRENTLY GENERATED. YOU'LL HAVE TO ADD THAT SEPARATELY 
#       AFTER COMBINING THE PLOTS.

# Species fixed differences
message("Plotting species differences...", stderr())
png(paste("fdiffs_files/", out_pref, "_", chr, "_fdiffs.png", sep = ""),
    width = 500, height = 3500)
  plot_spp_diverge(fdiff_tabname = spp_diffs_name, fst_tabname = fst_name,
                   chr = chr, chr_size = chr_size,
                   outfile_name = paste("fdiffs_files/", out_pref, "_", chr,
                                        "_fdiffs.bed", sep = ""))
dev.off()

# Recombination landscape
message("Plotting recombination...", stderr())
png(paste("recomb_files/", out_pref, "_", chr, "_recomb.png", sep = ""),
    width = 500, height = 3500)
  plot_recomb(recomb_tabname = recomb_name, ld_tabname = ld_name, chr = chr,
              chr_size = chr_size, 
              high_outname = paste("recomb_files/hotspots/", out_pref, "_", chr,
                                   "_recomb_hotspots.bed", sep = ""),
              low_outname = paste("recomb_files/suppressed/", out_pref, "_",
                                  chr, "_recomb_suppr.bed", sep = ""))
dev.off()

# Ancestry_HMM introgression
message("Plotting Ancestry_HMM results...", stderr())
png(paste("ahmm_files/", out_pref, "_", chr, "_ahmm.png", sep = ""),
    width = 500, height = 3500)
  plot_ahmm(ahmm_tabnames = ahmm_list,
            ahmm_outname = paste("ahmm_files/", out_pref, "_", chr, 
                                 "_ahmm.bed", sep = ""),
            dxy_tabname = dxy_name,
            scan_outname = paste("ahmm_files/gscan/", out_pref, "_", chr,
                                 "_gscan_intr.bed", sep = ""),
            chr = chr, chr_size = chr_size)
dev.off()

# ELAI introgression [WAITING FOR THIS TO FINISH]
message("Plotting ELAI results...", stderr())
png(paste("elai_files/", out_pref, "_", chr, "_elai.png", sep = ""),
    width = 500, height = 3500)
  plot_elai(elai_dose_file = elai_dose_name,
            elai_site_file = elai_snpinfo_name, elai_samples = elai_samples,
            elai_outname = paste("elai_files/", out_pref, "_", chr,
                                 "_elai.bed", sep = ""),
            dxy_tabname = dxy_name, chr = chr, chr_size = chr_size)
dev.off()

# Selection landscape
message("Plotting selection...", stderr())
png(paste("sel_files/", out_pref, "_", chr, "_sel.png", sep = ""),
    width = 500, height = 3500)
  plot_sel(tajima_tabname = tajd_name, ld_infile = ld_name,
           recomb_infile = recomb_name, chr = chr, chr_size = chr_size,
           bal_outname = paste("sel_files/balancing/", out_pref, "_", chr,
                               "_balance_sel.bed", sep = ""),
           dir_outname = paste("sel_files/direction/", out_pref, "_", chr,
                               "_direction_sel.bed", sep = ""),
           rec_ld_outname = paste("sel_files/direction/", out_pref, "_", chr,
                                  "_ld_rec_sel.bed", sep = ""))
dev.off()
