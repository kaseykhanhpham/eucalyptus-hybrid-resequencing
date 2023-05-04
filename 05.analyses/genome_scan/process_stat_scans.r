# R script to process egglib pi genome scan output
# Usage: Rscript process_pi_scans.r [infile] [stat] [sd_cutoff] [direction] [mean - OPTIONAL] [sd - OPTIONAL]
# infile: table output of egglib or vcftools to be processed for outlier windows
# stat: the stat of interest, should match the column reporting the stat for the window
# sd_cutoff: multiple of sd to set cutoff
# direction: set cutoff above or below mean
# mean: OPTIONAL - the mean to compare against
# sd: OPTIONAL - the sd to compare against. If you provide the mean you must provide this too.

# import command line arguments and read file
infile_name <- commandArgs(trailingOnly = TRUE)[1]
stat <- commandArgs(trailingOnly = TRUE)[2]
sd_cutoff <- commandArgs(trailingOnly = TRUE)[3]
direction <- commandArgs(trailingOnly = TRUE)[4]
in_table <- read.table(infile_name, header = TRUE, sep = "\t", na.strings="")

# Report general stats to stdout
print(paste("Sample", infile_name, ":"))
print(summary(in_table[, stat], na.rm = TRUE))
print(paste("sd:", SD = sd(in_table[, stat], na.rm = TRUE)))
# Print stat distribution to PNG file
png(paste(infile_name, "_", stat, "_distr.png", sep = ""), width = 500, height = 500)
    hist(in_table[, stat])
dev.off()

# Calculate summary stats
stat_mean <- mean(in_table[, stat], na.rm = TRUE)
stat_sd <- sd(in_table[, stat], na.rm = TRUE)

# Replace with user-provided mean and sd if given
if(length(commandArgs(trailingOnly = TRUE)) > 4){
    stat_mean <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
    stat_sd <- as.numeric(commandArgs(trailingOnly = TRUE)[6])
}

# Get outlier windows

flagged <- switch(direction,
    above = which(in_table[, stat] > (mean(stat_mean + 2*stat_sd))),
    below = which(in_table[, stat] < (mean(stat_mean - 2*stat_sd))))
flagged_windows <- in_table[flagged, c("chr", "start", "end")]

write.table(flagged_windows, paste(infile_name, ".flagged_", stat, ".tab", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)