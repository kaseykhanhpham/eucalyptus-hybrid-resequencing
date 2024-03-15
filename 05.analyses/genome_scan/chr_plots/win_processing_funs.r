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
    # get start and end indices and shift by 1 to match position indexing 
    # of other programs
    start_inds <- as.numeric(unlist(lapply(strsplit(inds, "-"), function(x) x[1]))) + 1
    end_inds <- as.numeric(unlist(lapply(strsplit(inds, "-"), function(x) x[2]))) + 1
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
    cutoff_point <- sort(infile[, "recombRate"])[round(nrow(infile) * cutoff)]
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
        cutoff_ci <- sort(infile[, "CI95LO"])[round(nrow(infile) * cutoff)]
        include_rows2 <- which(infile[, "CI95LO"] > cutoff_ci)
    } else if(dir == "below"){
        cutoff_ci <- sort(infile[, "CI95HI"])[round(nrow(infile) * cutoff)]
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

get_ahmm_windows <- function(infile_vec, chr, post_thresh, inds_thresh){
    # Purpose: get positions above a given posterior value for E. cordata ancestry (0,2)
    # Input: a vector of names of input files to process from AncestryHMM results, a
    #        minimum posterior to consider a site, and a minimum number of individuals in
    #        which the site is present 
    # Returns: a table in BED format

    # import input files
    infile_list <- lapply(infile_vec, function(infile_name) read.table(infile_name, sep = "\t", skip = 1, header = FALSE, col.names = c("chrom", "position", "hom_glob", "het", "hom_cord")))
    # subset to specified chromosome
    infile_list <- lapply(infile_list, function(infile) infile[which(infile$chrom == chr),])
    names(infile_list) <- infile_vec
    # get rows with a posterior larger than or equal to the specified threshold for (0,2)
    infile_filter <- lapply(infile_list, function(infile) which(infile$hom_cord > post_thresh))
    names(infile_filter) <- infile_vec
    
    # filter sites by how many individuals possess them
    final_sites_index <- c()
    pass_sites_inds <- c()
    for(i in c(1:nrow(infile_list[[1]]))){
        ind_site_passes <- unlist(lapply(infile_filter, function(infile) i %in% infile))
        # check that there are at least as many individuals with the site as the given threshold
        if(length(which(ind_site_passes)) >= inds_thresh){
            # add site to final site listing
            final_sites_index <- c(final_sites_index, i)
            # add inds to final ind list for each site
            pass_sites_inds <- c(pass_sites_inds, paste(infile_vec[which(ind_site_passes)], collapse = ","))
        }
    }
    # return BED of sites
    final_chr <- infile_list[[1]][final_sites_index, "chrom"]
    final_start <- infile_list[[1]][final_sites_index, "position"]
    final_end <- final_start + 1
    final_tab <- data.frame(chr = final_chr, start = final_start, end = final_end, inds = pass_sites_inds)
    return(final_tab)
}

get_elai_windows <- function(dose_name, sites_name, chr, dose_thresh, inds_thresh, sample_vec){
    # Purpose: get positions above a given posterior value for E. cordata ancestry dosage
    # Input: a vector of names of input files to process from Easy Local Ancestry Inference
    #        results, a minimum posterior to consider a site, a minimum number of 
    #        individuals in which the site is present, and a vector with the name of
    #        the samples
    # Returns: a table in BED format

    # import sites file
    sites_tab <- read.table(sites_name, header = TRUE, sep = "\t")
    # parse chromosomes
    site_chrs <- unlist(strsplit(sites_tab$rs, ":"), function(x) x[1])
    # import dosage file
    dose_tab <- scan(dose_name)
    dim(dose_tab) <- c(2,nrow(sites_tab),20) # 2 ancestral populations, num snps, 20 inds
    dose_tab <- dose_tab[2,,] # subset to cordata dosage only
    # subset by specified chromosome
    chr_filter <- which(sites_chr == chr)
    sites_tab <- sites_tab[chr_filter,]
    dose_tab <- dose_tab[chr_filter,]

    # get sites for each ind with dosage at least threshold
    dose_fil <- apply(dose_tab, 2, function(site) site > dose_thresh)
    
    # get sites for which dosage passed the threshold for the threshold number of individuals
    final_sites <- c()
    final_inds <- c()
    for(i in c(1:nrow(dose_tab))){
        if(length(which(dose_fil[i,])) >= inds_thresh){
            final_sites <- c(final_sites, i)
            # record individuals that passed the threshold filters
            final_inds <- paste(sample_vec[which(dose_fil[i])], sep = ",")
        }
    }
    # return BED format table of passed sites
    final_chrs <- rep(chr, length(final_sites))
    final_start <- sites_tab[final_sites, "pos"]
    final_end <- final_start + 1
    final_tab <- data.frame(chr = final_chrs, start = final_start, end = final_end, inds = final_inds)
    return(final_tab)
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
         density = -1, col = "#ea5a80", border = NA)
    # plot bars for protected regions
    rect(xleft = rep(4, nrow(protected_regions_merged)), ybottom = protected_regions_merged$start,
         xright = rep(6, nrow(protected_regions_merged)), ytop = protected_regions_merged$end,
         density = -1, col = "#501392", border = NA)
    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black")
}

plot_intr <- function(ahmm_filenames, ahmm_outname, elai_dose_file, elai_site_file, elai_samples, elai_outname, dxy_filename, dsuite_filename, chr, chr_size){
    # Plot introgression windows
    # Get windows of interest
    # AHMM
    ahmm_windows <- get_ahmm_windows(ahmm_filenames, chr, 0.95, 5)
    write.table(ahmm_windows, ahmm_outname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    ahmm_windows_merged <- merge_windows(ahmm_windows, 100)
    # ELAI
    elai_windows <- get_elai_windows(elai_dose_file, elai_site_file, chr, 1.75, 5, elai_samples)
    write.table(elai_windows, elai_outname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    elai_windows_merged <- merge_windows(elai_windows, 100)
    # dXY
    dxy_windows <- get_dxy_windows(dxy_filename, "glob_MR", "cord_MR", chr, 40, 0.50, "below")
    # fdm
    fdm_windows <- get_dsuite_windows(dsuite_filename, chr, "f_dM", 0.75, "above")

    # Get overlap between windows
    ahmm_intr_windows <- get_overlap(get_overlap(ahmm_windows_merged, dxy_windows, chr), fdm_windows, chr)
    elai_intr_windows <- get_overlap(get_overlap(elai_windows_merged, dxy_windows, chr), fdm_windows, chr)

    # merge windows for final plotting
    ahmm_intr_windows_merged <- merge_windows(ahmm_intr_windows, 1000)
    elai_intr_windows_merged <- merge_windows(elai_intr_windows, 1000)
    both_intr_windows <- get_overlap(ahmm_intr_windows_merged, elai_intr_windows_merged, chr)

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*chr_size*0.05), round(chr_size + chr_size*0.05)), 
         col = "white", xlab = "", ylab = "Position (bp)")
    # plot bars for AHMM regions
    rect(xleft = rep(4, nrow(ahmm_intr_windows_merged)), ybottom = ahmm_intr_windows_merged$start,
         xright = rep(6, nrow(ahmm_intr_windows_merged)), ytop = ahmm_intr_windows_merged$end,
         density = -1, col = "#ffb44a", border = NA)
    # plot bars for ELAI regions
    rect(xleft = rep(4, nrow(elai_intr_windows_merged)), ybottom = elai_intr_windows_merged$start,
         xright = rep(6, nrow(elai_intr_windows_merged)), ytop = elai_intr_windows_merged$end,
         density = -1, col = "#a85605", border = NA)
    # plot bars for overlaps
    rect(xleft = rep(4, nrow(both_intr_windows)), ybottom = both_intr_windows$start,
         xright = rep(6, nrow(both_intr_windows)), ytop = both_intr_windows$end,
         density = -1, col = "#000000", border = NA)
    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black")
}

plot_sel <- function(tajima_filename, ld_infile, recomb_infile, chr, chr_size){
    # Plot regions of selection (very high or low Tajima's D, increased LD, not low recombination rate)
    tajd_windows_high <- get_tajimad_windows(tajima_filename, chr, 40, 0.95, "above", 50000)
    tajd_windows_low <- get_tajimad_windows(tajima_filname, chr, 40, 0.05, "below", 50000)
    ld_windows <- get_ld_windows(ld_infile, chr, 0.75, "above")
    recomb_windows <- get_recomb_windows(recomb_infile, chr, 0.50, "above")

    # Overlap and merge windows
    recomb_windows_bal <- get_overlap(get_overlap(tajd_windows_high, ld_windows, chr), recomb_windows, chr)
    recomb_windows_bal_merge <- merge_windows(recomb_windows_bal, 1000)
    recomb_windows_dir <- get_overlap(get_overlap(tajd_windows_low, ld_windows, chr), recomb_windows, chr)
    recomb_windows_dir_merge <- merge_windows(recomb_windows_dir, 1000)

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*chr_size*0.05), round(chr_size + chr_size*0.05)), 
         col = "white", xlab = "", ylab = "Position (bp)")
    # plot bars for balancing selection
    rect(xleft = rep(4, nrow(recomb_windows_bal_merge)), ybottom = recomb_windows_bal_merge$start,
         xright = rep(6, nrow(recomb_windows_bal_merge)), ytop = recomb_windows_bal_merge$end,
         density = -1, col = "#1ddda3", border = NA)
    # plot bars for directional selection
    rect(xleft = rep(4, nrow(recomb_windows_dir_merge)), ybottom = recomb_windows_dir_merge$start,
         xright = rep(6, nrow(recomb_windows_dir_merge)), ytop = recomb_windows__merge$end,
         density = -1, col = "#008a60", border = NA)
    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black")
}
