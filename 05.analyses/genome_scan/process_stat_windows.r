### R script to process egglib stat genome scan output
### Usage: Rscript process_stat_windows.r [infile] [mode] [cols] [cutoff] [direction] [min_snps] [out]
###    infile: table output of genome scan program to be processed for outlier windows
###    mode: sd or percent, determines which measure cutoff is determined by
###    cols: the comma-separated indices of columns for chr,win_start,win_end,snp_count,statistic
###          where -1 indicates N/A value
###    sd_cutoff: multiple of sd / percentage at which to set cutoff
###    direction: set cutoff above or below mean / top or bottom percent
###    min_snps: The minimum number of SNPs to include a window in calculations
###    out: name of file to write with outlier windows

options(scipen = 999)

# import command line arguments and read file
infile_name <- commandArgs(trailingOnly = TRUE)[1]
mode <- commandArgs(trailingOnly = TRUE)[2]
col_nums <- commandArgs(trailingOnly = TRUE)[3]
cutoff <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
direction <- commandArgs(trailingOnly = TRUE)[5]
min_snps <- as.numeric(commandArgs(trailingOnly = TRUE)[6])
outfile_name <- commandArgs(trailingOnly = TRUE)[7]

in_table <- read.table(infile_name, header = TRUE, sep = "\t", na.strings = c("NA", "NaN", "nan"), as.is = TRUE)

chr_col <- as.numeric(strsplit(col_nums, ",")[[1]][1])
win_start_col <- as.numeric(strsplit(col_nums, ",")[[1]][2])
win_end_col <- as.numeric(strsplit(col_nums, ",")[[1]][3])
snp_count_col <- as.numeric(strsplit(col_nums, ",")[[1]][4])
stat_col <- as.numeric(strsplit(col_nums, ",")[[1]][5])

# convert relevant columns' classes
in_table[, win_start_col] <- as.numeric(in_table[, win_start_col])
if(win_end_col != -1){
    in_table[, win_end_col] <- as.numeric(in_table[, win_end_col])
}
if(snp_count_col != -1){
    in_table[, snp_count_col] <- as.numeric(in_table[, snp_count_col])
}
in_table[, stat_col] <- as.numeric(in_table[, stat_col])

# filter input by num SNPs
if(snp_count_col == -1){
    filtered_table <- in_table
} else {
    filtered_table <- in_table[which(in_table[, snp_count_col] >= min_snps),]
}

# Flag using standard deviation cutoff
if(mode == "sd"){
    # Calculate summary stats
    stat_mean <- mean(filtered_table[, stat_col], na.rm = TRUE)
    stat_sd <- sd(filtered_table[, stat_col], na.rm = TRUE)

    # Calculate cutoff
    sd_cutoff <- cutoff*stat_sd


    # Get outlier windows
    if(direction == "above"){
        flagged <- which(filtered_table[, stat_col] > (stat_mean + sd_cutoff))
    } else if(direction == "below"){
        flagged <- which(filtered_table[, stat_col] < (stat_mean - sd_cutoff))
    } else {
        warning(paste("Direction argument", direction, "not recognized. Specify above or below."))
    }

# Flag using percentage cutoff
} else if(mode == "percent"){
    # Calculate cutoff
    if(direction == "above"){
        cutoff_index <- round((1 - cutoff)*nrow(filtered_table))
    } else if(direction == "below"){
        cutoff_index <- round(cutoff*nrow(filtered_table))
    } else {
        warning(paste("Direction argument", direction, "not recognized. Specify above or below."))
    }
    percent_cutoff <- sort(filtered_table[, stat_col])[cutoff_index]

    # Get outlier windows
    if(direction == "above"){
        flagged <- which(filtered_table[, stat_col] >= percent_cutoff)
    } else if(direction == "below"){
        flagged <- which(filtered_table[, stat_col] <= percent_cutoff)
    }

} else {
    warning(paste("Mode value", mode, "not recognized. Specify sd or percent."))
}

# Build output table
if(win_end_col == -1){
    # If there isn't a window end column in the file, 
    # assume the window ends just before the next starts
    win_end <- unlist(sapply(c(1:(nrow(filtered_table) - 1)), function(i) filtered_table[i + 1, win_start_col] - 1))
    # assume the final end position is of the same step size as the penultimate one.
    last_step_size <- win_end[nrow(filtered_table) - 1] - filtered_table[nrow(filtered_table) - 1, win_start_col]
    last_end_pos <- filtered_table[nrow(filtered_table), win_start_col] + last_step_size
    win_end <- c(win_end, last_end_pos)

    flagged_windows <- filtered_table[flagged, c(chr_col, win_start_col)]
    flagged_windows <- cbind(flagged_windows, win_end[flagged])
    colnames(flagged_windows) <- c("Chr", "Start", "End")
} else {
    flagged_windows <- filtered_table[flagged, c(chr_col, win_start_col, win_end_col)]
    colnames(flagged_windows) <- c("Chr", "Start", "End")
}

write.table(flagged_windows, outfile_name, quote = FALSE, row.names = FALSE, col.names = TRUE)