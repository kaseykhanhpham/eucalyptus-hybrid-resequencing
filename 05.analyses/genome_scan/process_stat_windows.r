# R script to process egglib stat genome scan output
# Usage: Rscript process_stat_windows.r [infile] [stat] [mode] [cutoff] [direction] [min_snps] [mean - OPTIONAL] [sd - OPTIONAL]
# infile: table output of egglib or vcftools to be processed for outlier windows
# stat: the stat of interest, should match the column reporting the stat for the window
# mode: sd (standard deviation-based) or percent (top or bottom percentage of windows)
# sd_cutoff: multiple of sd / percentage at which to set cutoff
# direction: set cutoff above or below mean / top or bottom percent
# min_snps: The minimum number of SNPs to include a window in calculations
# mean: OPTIONAL - the mean to compare against
# sd: OPTIONAL - the sd to compare against. If you provide the mean you must provide this too.

# import command line arguments and read file
infile_name <- commandArgs(trailingOnly = TRUE)[1]
stat <- commandArgs(trailingOnly = TRUE)[2]
mode <- commandArgs(trailingOnly = TRUE)[3]
cutoff <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
direction <- commandArgs(trailingOnly = TRUE)[5]
min_snps <- as.numeric(commandArgs(trailingOnly = TRUE)[6])
in_table <- read.table(infile_name, header = TRUE, sep = "\t", na.strings="")

# filter input by num SNPs
filtered_table <- in_table[which(in_table$num_var >= min_snps),]

# Report general stats to stdout
write(paste("File:", infile_name, ":"), stdout())
write(summary(filtered_table[, stat]), stdout())
write(paste("sd:", SD = sd(filtered_table[, stat])), stdout())
# Print stat distribution to PNG file
png(paste(infile_name, "_distr.png", sep = ""), width = 500, height = 500)
    hist(filtered_table[, stat])
dev.off()

# Flag using standard deviation cutoff
if(mode == "sd"){
    # Calculate summary stats
    stat_mean <- mean(filtered_table[, stat])
    stat_sd <- sd(filtered_table[, stat])
    # Replace with user-provided mean and sd if given
    if(length(commandArgs(trailingOnly = TRUE)) > 6){
        stat_mean <- as.numeric(commandArgs(trailingOnly = TRUE)[7])
        stat_sd <- as.numeric(commandArgs(trailingOnly = TRUE)[8])
    }

    # Calculate cutoff
    sd_cutoff <- cutoff*stat_sd

    # Get outlier windows
    flagged <- switch(direction,
        above = which(filtered_table[, stat] > (stat_mean + sd_cutoff)),
        below = which(filtered_table[, stat] < (stat_mean - sd_cutoff)))

# Flag using percentage cutoff
} else if(mode == "percent"){
    # Calculate cutoff
    cutoff_index <- switch(direction,
        above = round((1 - cutoff)*nrow(filtered_table)),
        below = round(cutoff*nrow(filtered_table))
    )
    percent_cutoff <- sort(filtered_table[,stat])[cutoff_index]

    # Get outlier windows
    flagged <- switch(direction,
        above = which(filtered_table[, stat] >= percent_cutoff),
        below = which(filtered_table[, stat] <= percent_cutoff)
    )
} else {
    write("Mode value not recognized.", stdout())
}

flagged_windows <- filtered_table[flagged, c("chr", "start", "end")]

out_name <- paste(infile_name, ".", mode, cutoff, ".flagged.tab", sep = "")
write.table(flagged_windows, out_name, quote = FALSE, row.names = FALSE, col.names = TRUE)