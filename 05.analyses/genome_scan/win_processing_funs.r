# A set of R functions for plotting the results of genome scans to a given
# chromosome

## ------------------- ##
## WORKHORSE FUNCTIONS ##
## ------------------- ##
# Functions used for multiple graphs that identify outlier windows for various stats
# and return a table in a BED-like format for graphing

get_dsuite_windows <- function(filename, chr, stat, cutoff, dir){
    # Purpose: retrieve regions of interest from Dsuite investigate output
    # Input: the name of the Dsuite file, focal chromosome, 
    #        the statistic, cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # calculate cutoff value
    cutoff_num <- sort(infile[, stat])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, stat] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, stat] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # filter by chromosome
    chr_filter <- which(infile$chr == chr)
    # return output table
    outtab <- infile[intersect(include_rows, chr_filter) ,c("chr","windowStart","windowEnd")]
    colnames(outtab) <- c("chr", "start", "end") # standardize headings for downstream
    return(outtab)
}

get_pi_windows <- function(filename, spp, chr, min_sites, cutoff, dir){
    # Purpose: retrieve regions of interest from pixy pi output
    # Input: name of pixy pi file, species compared, focal chromosome,
    #        minimum sites to include a row, cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # filter by species and min sites
    infile <- infile[which(infile$no_sites >= min_sites),]
    infile <- infile[which(infile$pop == spp),]
    # calculate cutoff value
    cutoff_num <- sort(infile[, "avg_pi"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_pi"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_pi"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # create chromosome filter
    chr_fil <- which(infile$chromosome == chr)
    # return output table
    outtab <- infile[intersect(include_rows, chr_filter), c("chromosome", "window_pos_1", "window_pos_2")]
    colnames(outtab) <- c("chr", "start", "end")
    return(outtab)
}

get_dxy_windows <- function(filename, spp1, spp2, chr, min_sites, cutoff, dir){
    # Purpose: retrieve regions of interest from pixy dxy output
    # Input: name of pixy dxy file, species compared, focal chromosome,
    #        minimum sites to include a row, cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # filter by species and min sites
    infile <- infile[which(infile$no_sites >= min_sites),]
    spp_rows <- c(which(infile$pop1 == spp1 & infile$pop2 == spp2), which(infile$pop1 == spp2 & infile$pop2 == spp1))
    infile <- infile[spp_rows,]
    # calculate cutoff value
    cutoff_num <- sort(infile[, "avg_dxy"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_dxy"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_dxy"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # create chromosome filter
    chr_filter <- which(infile$chromosome == chr)
    # return output table
    outtab <- infile[intersect(include_rows, chr_filter), c("chromosome", "window_pos_1", "window_pos_2")]
    colnames(outtab) <- c("chr", "start", "end")
    return(outtab)
}

get_fst_windows <- function(filename, spp1, spp2, chr, min_sites, cutoff, dir){
    # Purpose: retrieve regions of interest from pixy dxy output
    # Input: name of pixy dxy file, species compared, focal chromosome,
    #        minimum sites to include a row, cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # filter by species and min sites
    infile <- infile[which(infile$no_snps >= min_sites),]
    spp_rows <- c(which(infile$pop1 == spp1 & infile$pop2 == spp2), which(infile$pop1 == spp2 & infile$pop2 == spp1))
    infile <- infile[spp_rows,]
    # calculate cutoff value
    cutoff_num <- sort(infile[, "avg_wc_fst"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_wc_fst"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_wc_fst"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # create chromosome filter
    chr_list <- which(infile$chromosome == chr)
    # return output table
    outtab <- infile[intersect(include_rows, chr_list), c("chromosome", "window_pos_1", "window_pos_2")]
    colnames(outtab) <- c("chr", "start", "end")
    return(outtab)
}

get_tajimad_windows <- function(filename, chr, min_sites, cutoff, dir, bin_size){
    # Purpose: retrieve regions of interest from vcftools Tajima's D output
    # Input: name of vcftools Tajima's D file, focal chromosome, 
    #        minimum sites to include a row, cutoff percentile,  direction,
    #        and size of window
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # filter by min sites
    infile <- infile[which(infile$N_SNPS >= min_sites),]
    # calculate cutoff value
    cutoff_num <- sort(infile[, "TajimaD"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "TajimaD"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "TajimaD"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # create chromosome and final filter
    chr_filter <- which(infile$CHROM == chr)
    final_rows <- intersect(chr_filter, include_rows)
    # shift indexing on positions to start at 1, such that it matches 
    # with other outputs
    start_inds <- infile[final_rows, "BIN_START"] + 1
    end_inds <- infile[final_rows, "BIN_START"] + bin_size
    # return output table
    outtab <- data.frame(chr = infile[final_rows, "CHROM"], start = start_inds, end = end_inds)
    return(outtab)
}

get_ld_windows <- function(filename, chr, cutoff, dir){
    # Purpose: retrieve regions of interest for LD from previous summary scripts
    # Input: name of LD window file from fit_rsq_curve.r, chromosome name,
    #        cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # parse start and end of windows from file_name column
    chrs <- unlist(lapply(strsplit(infile$file_name, "_"), function(x) x[1]))
    inds <- unlist(lapply(strsplit(infile$file_name, "_"), function(x) x[2]))
    start_inds <- unlist(lapply(strsplit(inds, "-"), function(x) x[1]))
    end_inds <- unlist(lapply(strsplit(inds, "-"), function(x) x[2]))
    # calculate cutoff value
    cutoff_num <- sort(infile[, "ld"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "ld"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "ld"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # create chromosome and final filter
    chr_filter <- which(chrs == chr)
    final_filter <- intersect(include_rows, chr_filter)
    # return output table
    outtab <- data.frame(chr = chrs[final_filter], start = start_inds[final_filter], end = end_inds[final_filter])
    outtab <- outtab[order(outtab$start),]
    return(outtab)
}

get_recomb_windows <- function(filename, chr, cutoff, dir){
    # Purpose: retrieve regions of interest for recombination from ReLERNN BSCORRECT
    # Input: name of ReLERNN window file from BSCORRECT output, focal chromosome, 
    #        cutoff percentile (applies to both point estimate and CI boundary), 
    #        and direction (applies to both point estimate and CI boundary)
    # Note: Assumes that a suitable minimum SNPs per window has already been set
    #       during ReLERNN run and that ReLERNN was run separately for each chr
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # calculate cutoff value for point est
    cutoff_point <- order(infile[, "recombRate"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows1 <- which(infile[, "recombRate"] > cutoff_point)
    } else if(dir == "below"){
        include_rows1 <- which(infile[, "recombRate"] < cutoff_point)
    } else {
        print("Please specify above or below for direction.")
        include_rows1 <- c()
    }
    # calculate cutoff value for confidence interval
    if(dir == "above"){
        cutoff_ci <- order(infile[, "CI95LO"])[round(nrow(infile) * cutoff)]
        include_rows2 <- which(infile[, "CI95LO"] > cutoff_ci)
    } else if(dir == "below"){
        cutoff_ci <- order(infile[, "CI95HI"])[round(nrow(infile) * cutoff)]
        include_rows2 <- which(infile[, "CI95HI"] < cutoff_ci)
    } else {
        include_rows2 <- c()
    }
    # create chromosome filter and final filter
    chr_filter <- which(infile$chrom == chr)
    final_filter <- intersect(chr_filter, intersect(include_rows1, include_rows2))
    # return output table
    outtab <- infile[final_filter, c("chrom", "start", "end")]
    colnames(outtab) <- c("chr", "start", "end")
    return(outtab)
}

merge_windows <- function(wintab, thresh){
    # Purpose: merge adjacent windows into a single window entry for a BED table
    # Input: a dataframe in BED format, the maximum amount of base pairs apart
    #        that two windows can be to merge them
    # Returns: a table in BED format

    # initialize cols for output dataframe
    chr_new <- c()
    start_new <- c()
    end_new <- c()

    # add dummy entry at bottom of input table to evaluate last real row in coming loop
    new_row <- c("DUMMY_ROW", 1, 2)
    wintab[(nrow(wintab) + 1),] <- new_row

    # initialize storage variables for loop
    curr_chr <- wintab[1, "chr"]
    curr_start <- wintab[1, "start"]
    curr_end <- wintab[1, "end"]

    # loop through input dataframe
    for(i in c(2:nrow(wintab))){
        switch_chr <- wintab[i, "chr"] != curr_chr
        bigger_thr <- wintab[i, "start"] - curr_end > thresh
        if(switch_chr | bigger_thr){
            # if switching chromosomes or regions are farther than threshold,
            # write and update all variables
            chr_new <- c(chr_new, curr_chr)
            start_new <- c(start_new, curr_start)
            end_new <- c(end_new, curr_end)
            
            curr_chr <- wintab[i, "chr"]
            curr_start <- wintab[i, "start"]
            curr_end <- wintab[i, "end"]
        } else {
            # if regions are closer than threshold, update endpoint 
            curr_end <- wintab[i, "end"]
        }
    }

    # return output dataframe
    outtab <- data.frame(chr = chr_new, start = start_new, end = end_new)
    return(outtab)
}

get_overlap <- function(tab1, tab2, chr){
    # get overlapping regions between two BED-like tables
    # expand input tables into individual sites
    tab1_sites <- unlist(apply(tab1, 1, function(x) c(x["start"]:x["end"])))
    tab2_sites <- unlist(apply(tab2, 1, function(x) c(x["start"]:x["end"])))

    # get intersection of sites between tabs
    int_sites <- intersect(tab1_sites, tab2_sites)

    # convert back into BED-like format
    # get start positions (positions where the previous entry is not 1 less)
    int_starts <- unlist(sapply(c(2:length(int_sites)), function(i) int_sites[i] != (int_sites[(i - 1)] + 1)))
    int_starts <- c(int_sites[1], int_starts)
    # get end positions (positions where the next entry is not 1 more)
    int_ends <- unlist(sapply(c(1:(length(int_sites) - 1)), function(i) int_sites[i] != (int_sites[(i + 1)] - 1)))
    int_ends <- c(int_ends, int_sites[length(int_sites)])
    # make chromosome column
    chrs <- rep(chr, length(int_starts))
    # construct output table
    outtab <- data.frame(chr = chrs, start = int_starts, end = int_ends)
    return(outtab)
}

## --------------------------- ##
## SPECIFIC GRAPHING FUNCTIONS ##
## --------------------------- ##
# Functions that graph blocks across chromosomes for a certain pattern
# (e.g., putative introgression, putative selective sweep)

plot_spp_diverge <- function(fdiff_tabname, fst_tabname, chr, chr_size){
    # Plot chromosome differences between E. globulus (reference) and E. cordata
    # on a given chromosome

    # get windows to plot
    fdiff_tab <- read.table(fdiff_tabname, header = TRUE, sep = "\t")
    fst_80_wins <- get_fst_windows(fst_tabname, "glob_ref", "cord_MR", chr, 15, 0.8, "above")
    fst_95_wins <- get_fst_windows(fst_tabname, "glob_ref", "cord_MR", chr, 15, 0.95, "above")
    # (lower minimum site requirement because FST only considers variable sites)

    # merge windows that are close to each other (< 1000 bp)
    fst_80_wins_merged <- merge_windows(fst_80_wins, 1000)
    fst_95_wins_merged <- merge_windows(fst_95_wins, 1000)

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*chr_size*0.05), round(chr_size + chr_size*0.05)), 
         col = "white", xlab = "", ylab = "Position (bp)")
    # Plot bars for FST > 0.80
    rect(xleft = rep(4, nrow(fst_80_wins_merged)), ybottom = fst_80_wins_merged$start, 
         xright = rep(6, nrow(fst_80_wins_merged)), ytop = fst_80_wins_merged$end,
         density = -1, col = "#ababab", border = NA)
    # Plot bars for FST > 0.95
    rect(xleft = rep(4, nrow(fst_95_wins_merged)), ybottom = fst_95_wins_merged$start, 
         xright = rep(6, nrow(fst_95_wins_merged)), ytop = fst_95_wins_merged$end,
         density = -1, col = "#494949", border = NA)
    # Plot fixed differences
    rect(xleft = rep(4, nrow(fdiff_tab)), ybottom = fdiff_tab$pos, 
         xright = rep(6, nrow(fdiff_tab)), ytop = fdiff_tab$pos,
         density = -1, col = "#000000", border = NA)
    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black")
}

plot_recomb <- function(recomb_tabname, ld_tabname, chr, chr_size){
    # Plot recombination hotspots and regions of recombination suppression
    # on a given chromosome
    
    # get hotspots and merge close windows (< 1kbp)
    recomb_high_wins <- get_recomb_windows(recomb_tabname, chr, 0.95, "above")
    recomb_high_wins_merged <- merge_windows(recomb_high_wins, 1000)
    # get protected regions and merge close windows (< 1 kbp)
    recomb_low_wins <- get_recomb_windows(recomb_tabname, chr, 0.05, "below")
    ld_high_wins <- get_ld_windows(ld_tabname, chr, 0.75, "above")
    protected_regions <- get_overlap(recomb_low_wins, ld_high_wins, chr)
    protected_regions_merged <- merge_windows(protected_regions, 1000)

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*chr_size*0.05), round(chr_size + chr_size*0.05)), 
         col = "white", xlab = "", ylab = "Position (bp)")
    # plot bars for hotspots
    rect(xleft = rep(4, nrow(recomb_high_wins_merged)), ybottom = recomb_high_wins_merged$start,
         xright = rep(6, nrow(recomb_high_wins_merged)), ytop = recomb_high_wins_merged$end,
         density = -1, col = "#f3981a", border = NA)
    # plot bars for protected regions
    rect(xleft = rep(4, nrow(protected_regions_merged)), ybottom = protected_regions_merged$start,
         xright = rep(6, nrow(protected_regions_merged)), ytop = protected_regions_merged$end,
         density = -1, col = "#646464", border = NA)
    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black")
}