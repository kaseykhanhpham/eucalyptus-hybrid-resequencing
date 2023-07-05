# R script to calculate average pi/dxy from a pixy output file given
# a corresponding BED file of windows (must match windows in pixy file)
#
# Usage: Rscript calc_window_stats.r [pixy file] [bed file] [min cov] [stat] [pop 1] [pop 2]
#        pixy file: the file name of the pixy file to calculate stats from
#        bed file: the bed file of windows from which to calculate stat
#        min cov: minimum number of sites to include window in overall calculation
#        stat: the stat to be calculated, either pi, dxy, or fst
#        pop 1: a comma-delimited list of populations to consider for calculation of stats
#        pop 2: only for dxy and fst. a comma-delimited list of populations paired to pop 1
#               to consider for calculation of stats

# get cli arguments
pixy_file_name <- commandArgs(trailingOnly = TRUE)[1]
bed_file_name <- commandArgs(trailingOnly = TRUE)[2]
min_cov <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
stat <- commandArgs(trailingOnly = TRUE)[4]
pop1 <- unlist(strsplit(commandArgs(trailingOnly = TRUE)[5], ","))
if(stat %in% c("dxy", "fst")){
    pop2 <- unlist(strsplit(commandArgs(trailingOnly = TRUE)[6], ","))
}

# import input files
pixy_file <- read.table(pixy_file_name, header = TRUE)
bed_file <- read.table(bed_file_name, header = FALSE)

# correct for index-1 in bedfile
bed_file[,2] <- bed_file[,2] - 1

# subset by bedfile-specified windows
pixy_bed_mask <- unlist(apply(bed_file, 1, function(row) which(pixy_file$chromosome == row[1] &
                                                               pixy_file$window_pos_1 == as.numeric(row[2]) &
                                                               pixy_file$window_pos_2 == as.numeric(row[3]))))
# subset by minimum coverage
if(stat %in% c("pi", "dxy")){
    pixy_cov_mask <- which(pixy_file$no_sites >= min_cov)
} else {
    pixy_cov_mask <- which(pixy_file$no_snps >= min_cov)
}
# subset by population
if(stat == "pi"){
    pixy_pop_mask <- which(pixy_file$pop %in% pop1)
} else {
    pixy_pop_mask <- unlist(sapply(seq(1,length(pop1)), function(i) which(pixy_file$pop1 == pop1[i] &
                                                                          pixy_file$pop2 == pop2[i])))
}
pixy_full_mask <- intersect(pixy_bed_mask, pixy_cov_mask, pixy_pop_mask)

# calculate mean stat over windows subset
if(stat %in% c("pi", "dxy")){
    stat_calc <- sum(pixy_file[pixy_full_mask,"count_diffs"]) /
                 sum(pixy_file[pixy_full_mask,"count_comparisons"])
} else {
    stat_calc <- mean(pixy_file[pixy_full_mask, "avg_wc_fst"], na.rm = TRUE)
}

# calculate SD forstats in window subset
if(stat %in% c("pi", "dxy")){
    stat_col <- paste("avg", stat, sep = "_")
} else {
    stat_col <- paste("avg", "wc", stat, sep = "_")
}
stat_sd <- sd(pixy_file[pixy_full_mask, stat_col], na.rm = TRUE)

# print histogram of stat in window subset
hist_name <- paste(pixy_file_name, bed_file, "hist.png", sep = "_")
png(hist_name, height = 1000, width = 1000, units = "px")
    hist(pixy_file[pixy_full_mask, stat_col])
dev.off()

# print to stdout()
write(paste("Stat Calculations for", pixy_file_name, "and", bed_file_name), "", append = TRUE)
write(paste("mean", stat, ":", stat_calc, "\n"), "", append = TRUE)
write(paste("sd of", stat, ":", stat_sd, "\n"), "", append = TRUE)
write(paste("hist written to:", hist_name, "\n"), "", append = TRUE)
