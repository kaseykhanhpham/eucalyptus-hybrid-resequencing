# A set of R functions for plotting the results of genome scans to a given
# chromosome

get_dsuite_windows <- function(filename, chr, stat, cutoff, dir){
    # Purpose: retrieve regions of interest from Dsuite investigate output
    # Input: the name of the Dsuite file, focal chromosome, 
    #        the statistic, cutoff percentile, and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # filter by chromosome
    infile <- infile[which(infile$chr == chr),]
    # calculate cutoff value
    cutoff_num <- order(infile[, stat])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, stat] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, stat] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # return output table
    outtab <- infile[include_rows,c("chr","windowStart","windowEnd")]
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
    # filter by chromosome, species, and min sites
    infile <- infile[which(infile$chromosome == chr),]
    infile <- infile[which(infile$no_sites >= min_sites),]
    infile <- infile[which(infile$pop == spp),]
    # calculate cutoff value
    cutoff_num <- order(infile[, "avg_pi"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_pi"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_pi"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # return output table
    outtab <- infile[include_rows, c("chromosome", "window_pos_1", "window_pos_2")]
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
    # filter by chromosome, species, and min sites
    infile <- infile[which(infile$chromosome == chr),]
    infile <- infile[which(infile$no_sites >= min_sites),]
    spp_rows <- c(which(infile$pop1 == spp1 & infile$pop2 == spp2), which(infile$pop1 == spp2 & infile$pop2 == spp1))
    infile <- infile[spp_rows,]
    # calculate cutoff value
    cutoff_num <- order(infile[, "avg_dxy"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_dxy"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_dxy"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # return output table
    outtab <- infile[include_rows, c("chromosome", "window_pos_1", "window_pos_2")]
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
    # filter by chromosome, species, and min sites
    infile <- infile[which(infile$chromosome == chr),]
    infile <- infile[which(infile$no_snps >= min_sites),]
    spp_rows <- c(which(infile$pop1 == spp1 & infile$pop2 == spp2), which(infile$pop1 == spp2 & infile$pop2 == spp1))
    infile <- infile[spp_rows,]
    # calculate cutoff value
    cutoff_num <- order(infile[, "avg_wc_fst"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "avg_wc_fst"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "avg_wc_fst"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # return output table
    outtab <- infile[include_rows, c("chromosome", "window_pos_1", "window_pos_2")]
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
    # filter by chromosome, species, and min sites
    infile <- infile[which(infile$CHROM == chr),]
    infile <- infile[which(infile$N_SNPS >= min_sites),]
    # calculate cutoff value
    cutoff_num <- order(infile[, "TajimaD"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "TajimaD"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "TajimaD"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # shift indexing on positions to start at 1, such that it matches 
    # with other outputs
    start_inds <- infile[include_row, "BIN_START"] + 1
    end_inds <- infile[include_row, "BIN_START"] + bin_size
    # return output table
    outtab <- data.frame(chr = infile[include_rows, "CHROM"], start = start_inds, end = end_inds)
    return(outtab)
}

get_ld_windows <- function(filename, cutoff, dir){
    # Purpose: retrieve regions of interest for LD from previous summary scripts
    # Input: name of LD window file from fit_rsq_curve.r, cutoff percentile,
    #        and direction
    # Returns: a table with equivalent structure to a simple BED file with
    #          the regions of interest

    infile <- read.table(filename, sep = "\t", header = TRUE)
    # parse start and end of windows from file_name column
    chrs <- unlist(lapply(strsplit(infile$file_name, "_"), function(x) x[1]))
    inds <- unlist(lapply(strsplit(infile$file_name, "_"), function(x) x[2]))
    start_inds <- unlist(lapply(strsplit(inds, "-"), function(x) x[1]))
    end_inds <- unlist(lapply(strsplit(inds, "-"), function(x) x[2]))
    # calculate cutoff value
    cutoff_num <- order(infile[, "ld"])[round(nrow(infile) * cutoff)]
    # get rows to include in output based on direction specified
    if(dir == "above"){
        include_rows <- which(infile[, "ld"] > cutoff_num)
    } else if(dir == "below"){
        include_rows <- which(infile[, "ld"] < cutoff_num)
    } else {
        print("Please specify above or below for direction.")
        include_rows <- c()
    }
    # return output table
    outtab <- data.frame(chr = chrs[include_rows], start = start_inds[include_rows], end = end_inds[include_rows])
    outtab <- outtab[sort(outtab$start),] # CHECK THIS
    return(outtab)
}

get_recomb_windows <- function(filename, cutoff, dir){
    # Purpose: retrieve regions of interest for recombination from ReLERNN BSCORRECT
    # Input: name of ReLERNN window file from BSCORRECT output, cutoff percentile
    #        (applies to both point estimate and CI boundary), and direction (applies
    #        to both point estimate and CI boundary)
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
    # return output table
    outtab <- infile[intersect(include_rows1, include_rows2), c("chrom", "start", "end")]
    colnames(outtab) <- c("chr", "start", "end")
    return(outtab)
}

merge_windows <- function(wintab, thresh){
    # Purpose: merge adjacent windows into a single window entry for a BED table
    # Input: a dataframe in BED format, the maximum amount of base pairs apart
    #        that two windows can be to merge them
    # Returns: a table in BED format
}