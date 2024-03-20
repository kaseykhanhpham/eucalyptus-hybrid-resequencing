# Plot output of `bcftools query` for VCF annotations
# Usage: Rscript plot_snp_stats.r -d [depth_file] -g [gq_file] -s [sp_file] -i [info_file] -o [out_prefix]
# Note: assumes all outputs of BCFTools query have headers, which are not generated by default by the program.

library(argparser)
parser <- arg_parser("A script to plot sequencing coverage")

parser <- add_argument(parser, "-d", help = "Depth output file from bcftools query", default = "query_depth.txt", type = "character")
parser <- add_argument(parser, "-g", help = "Genotype Quality output file from bcftools query", default = "query_gq.txt", type = "character")
parser <- add_argument(parser, "-s", help = "Strand P-value output file from bcftools query", default = "query_sp.txt", type = "character")
parser <- add_argument(parser, "-i", help = "Info output file from bcftools query", default = "query_info.txt", type = "character")
parser <- add_argument(parser, "-o", help = "Output prefix string", default = "out", type = "character")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

dp <- read.table(args_list[["d"]], header = TRUE, na.strings=".")
total_dp <- apply(dp, 1, sum)
rm(dp)

png(filename = paste(args_list[["o"]], "dp_distr.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    hist(total_dp, breaks = 50, main = paste(args_list[["o"]], "total depth distribution"), xlab = "Total DP")
dev.off()

gq <- read.table(args_list[["g"]], header = TRUE, na.strings=".")
avg_gq <- apply(gq, 1, mean)
rm(gq)

png(filename = paste(args_list[["o"]], "gq_distr.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    hist(avg_gq, breaks = 50, main = paste(args_list[["o"]], "average genotype quality distribution"), xlab = "Average GQ")
dev.off()

png(filename = paste(args_list[["o"]], "dp_vs_gq.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    plot(x = avg_gq, y = total_dp, main = paste(args_list[["o"]], "depth/genotype quality correlation"), xlab = "GQ", ylab = "DP")
dev.off()

rm(avg_gq)

sp <- read.table(args_list[["s"]], header = TRUE, na.strings=".")
avg_sp <- apply(sp, 1, mean)
rm(sp)

png(filename = paste(args_list[["o"]], "sp_distr.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    hist(avg_sp, breaks = 50, main = paste(args_list[["o"]], "average strand bias distribution"), xlab = "Average SP")
dev.off() 

rm(avg_sp)

info <- read.table(args_list[["i"]], header = TRUE, na.strings=".", colClasses = c("character", "numeric", "numeric", "numeric", "numeric", "numeric"))

png(filename = paste(args_list[["o"]], "qual_distr.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    hist(info$QUAL, breaks = 50, main = paste(args_list[["o"]], "variant quality distribution"), xlab = "QUAL")
dev.off()

png(filename = paste(args_list[["o"]], "dp_vs_qual.png", sep = "_"),
    width = 1000, height = 1000, unit = "px")
    plot(x = info$QUAL, y = total_dp, main = paste(args_list[["o"]], "depth/quality correlation"), xlab = "QUAL", ylab = "DP")
dev.off()

write(paste("done plotting", args_list[["o"]]), stdout())