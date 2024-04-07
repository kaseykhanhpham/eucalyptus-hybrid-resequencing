# Functions to plot Local Ancestry Inference with lower threshold for sharing between individuals
# than in win_processing_funs.r

plot_ahmm_low <- function(ahmm_tabnames, ahmm_outname, dxy_tabname, dsuite_tabname, scan_outname, chr, chr_size){
    # Plot introgression windows for genome scan statistics and Ancestry_HMM
    max_chr_size <- 70214608
    # Get windows of interest
    # AHMM
    ahmm_seeds <- get_ahmm_windows(ahmm_tabnames, chr, 0.95, 4)
    if(nrow(ahmm_seeds) > 0){ # only output and merge table if there are sites to output
        write.table(ahmm_seeds, paste("seeds", ahmm_outname, sep = "_"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        ahmm_seeds_merged <- merge_windows(ahmm_seeds[,c(1:3)], 500, chr)
    } else {
        ahmm_seeds_merged <- data.frame()
    }
    # Expand AHMM regions
    # get areas where posterior > 0.89
    ahmm_surr <- get_ahmm_windows(ahmm_tabnames, chr, 0.89, 4)
    if(nrow(ahmm_surr) > 0){
        ahmm_surr_merged <- merge_windows(ahmm_surr[,c(1:3)], 500, chr)
    } else { 
        ahmm_surr_merged <- data.frame()
    }
    # subset to areas where posterior > 0.89 AND contains at least one region where posterior > 0.95
    final_start <- c()
    final_end <- c()
    if(nrow(ahmm_seeds_merged) > 0){
        for(i in c(1:nrow(ahmm_seeds_merged))){
            lower <- ahmm_surr_merged$start <= ahmm_seeds_merged[i, "start"]
            higher <- ahmm_surr_merged$end >= ahmm_seeds_merged[i, "end"]
            if(any(lower & higher)){ # check if there are any regions with posterior > 0.89 that contain
                                      # a seed
                # if there are, add that region to the final set
                final_start <- c(final_start, ahmm_surr_merged[which(lower & higher), "start"])
                final_end <- c(final_end, ahmm_surr_merged[which(lower & higher), "end"])
            } else {
                # if not, add original seed location
                final_start <- c(final_start, ahmm_seeds_merged[i, "start"])
                final_end <- c(final_end, ahmm_seeds_merged[i, "end"])
            }
        }
    }
    # construct final AHMM table of 0.95 posterior regions and surrounding 0.89 posterior regions
    ahmm_windows <- data.frame(chr = rep(chr, length(final_start)), start = final_start, end = final_end)
    ahmm_windows_merged <- merge_windows(ahmm_windows, 5001, chr)
    if(nrow(ahmm_windows) > 0){
        write.table(ahmm_windows_merged, paste("merged", ahmm_outname, sep = "_"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }

    # dXY
    dxy_windows <- get_dxy_windows(dxy_tabname, "glob_MR", "cord_MR", chr, 40, 0.10, "below")
    dxy_windows_merged <- merge_windows(dxy_windows, 5001, chr)
    # fdM
    fdm_windows <- get_dsuite_windows(dsuite_tabname, chr, "f_dM", 0.90, "above")
    fdm_windows_merged <- merge_windows(fdm_windows, 5001, chr)

    # Get overlap between genome scan windows
    dxy_fdm_overl <- get_overlap(dxy_windows_merged, fdm_windows_merged, chr)
    if(nrow(dxy_fdm_overl) > 0){ # plot overlap
        write.table(dxy_fdm_overl, scan_outname, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*max_chr_size*0.05), round(max_chr_size + max_chr_size*0.05)), 
         col = "white", xlab = "", xaxt = "n", ylab = "Position (bp)")
    # plot bars for dxy regions
    if(nrow(dxy_windows_merged) > 0){
        rect(xleft = rep(4, nrow(dxy_windows_merged)), ybottom = dxy_windows_merged$start,
            xright = rep(6, nrow(dxy_windows_merged)), ytop = dxy_windows_merged$end,
            density = -1, col = "#18CC3E", border = "#18CC3E", lwd = 1.5)
    }
    # plot bars for fdm regions
    if(nrow(fdm_windows_merged) > 0){ # check that table is populated for this chromosome
            rect(xleft = rep(4, nrow(fdm_windows_merged)), ybottom = fdm_windows_merged$start,
                xright = rep(6, nrow(fdm_windows_merged)), ytop = fdm_windows_merged$end,
                density = -1, col = "#FFF147", border = "#FFF147", lwd = 1.5)
    }
    # plot bars for overlaps
    if(nrow(dxy_fdm_overl) > 0){
        rect(xleft = rep(4, nrow(dxy_fdm_overl)), ybottom = dxy_fdm_overl$start,
            xright = rep(6, nrow(dxy_fdm_overl)), ytop = dxy_fdm_overl$end,
            density = -1, col = "#b5d42a", border = "#b5d42a", lwd = 1.5)
    }
    # plot bars for AHMM regions
    if(nrow(ahmm_windows_merged) > 0){ # check that table is populated for this chromosome
        rect(xleft = rep(4, nrow(ahmm_windows_merged)), ybottom = ahmm_windows_merged$start,
            xright = rep(6, nrow(ahmm_windows_merged)), ytop = ahmm_windows_merged$end,
            density = -1, col = "#263CA9", border = "#263CA9", lwd = 1.5)
    }

    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black", lwd = 4)
}

# I would not worry too much about ELAI, there might not even be any sites to graph.
plot_elai_low <- function(elai_dose_file, elai_site_file, elai_samples, elai_outname, dxy_tabname, dsuite_tabname, chr, chr_size){
    # Plot introgression windows
    max_chr_size <- 70214608
    # Get windows of interest
    # ELAI
    elai_seeds <- get_elai_windows(elai_dose_file, elai_site_file, chr, 1.75, 4, elai_samples)
    if(nrow(elai_seeds) > 0){ # only output table if there are sites to output
        write.table(elai_seeds, paste("seed", elai_outname, sep = "_"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        elai_seeds_merged <- merge_windows(elai_seeds[,c(1:3)], 500, chr)
    } else {
        elai_seeds_merged <- data.frame()
    }
    # Expand ELAI regions
    # get areas where dose > 1.5
    elai_surr <- get_elai_windows(elai_dose_file, elai_site_file, chr, 1.5, 4, elai_samples)
    if(nrow(elai_surr) > 0){
        elai_surr_merged <- merge_windows(elai_surr[,c(1:3)], 500, chr)
    } else { 
        elai_surr_merged <- data.frame()
    }
    # subset to areas where dosage > 1.5 AND contains at least one region where dosage > 1.75
    final_start <- c()
    final_end <- c()
    if(nrow(elai_seeds_merged) > 0){
        for(i in c(1:nrow(elai_seeds_merged))){
            lower <- elai_surr_merged$start <= elai_seeds_merged[i, "start"]
            higher <- elai_surr_merged$end >= elai_seeds_merged[i, "end"]
            if(any(lower & higher)){ # check if there are any regions with dosage > 1.5 that contain
                                      # a seed
                # if there are, add that region to the final set
                final_start <- c(final_start, elai_surr_merged[which(lower & higher), "start"])
                final_end <- c(final_end, elai_surr_merged[which(lower & higher), "end"])
            } else {
                # if not, add original seed location
                final_start <- c(final_start, elai_seeds_merged[i, "start"])
                final_end <- c(final_end, elai_seeds_merged[i, "end"])
            }
        }
    }
    # construct final ELAI table of 1.75 dosage regions and surrounding 1.5 dosage regions
    elai_windows <- data.frame(chr = rep(chr, length(final_start)), start = final_start, end = final_end)
    elai_windows_merged <- merge_windows(elai_windows, 5001, chr)
    if(nrow(elai_windows) > 0){
        write.table(elai_windows_merged, paste("merged", elai_outname, sep = "_"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }

    # dXY
    dxy_windows <- get_dxy_windows(dxy_tabname, "glob_MR", "cord_MR", chr, 40, 0.10, "below")
    dxy_windows_merged <- merge_windows(dxy_windows, 5001, chr)
    # fdM
    fdm_windows <- get_dsuite_windows(dsuite_tabname, chr, "f_dM", 0.90, "above")
    fdm_windows_merged <- merge_windows(fdm_windows, 5001, chr)

    # Get overlap between genome scan windows
    dxy_fdm_overl <- get_overlap(dxy_windows_merged, fdm_windows_merged, chr)
    # don't need to write output because already done in plot_ahmm

    # make canvas of the correct size (adding 5% padding to either end of the y-axis)
    plot(x = c(0, 10), y = c(round(-1*max_chr_size*0.05), round(max_chr_size + max_chr_size*0.05)), 
         col = "white", xlab = "", xaxt = "n", ylab = "Position (bp)")
    # plot bars for dxy regions
    if(nrow(dxy_windows_merged) > 0){
        rect(xleft = rep(4, nrow(dxy_windows_merged)), ybottom = dxy_windows_merged$start,
            xright = rep(6, nrow(dxy_windows_merged)), ytop = dxy_windows_merged$end,
            density = -1, col = "#18CC3E", border = "#18CC3E", lwd = 1.5)
    }
    # plot bars for fdm regions
    if(nrow(fdm_windows_merged) > 0){ # check that table is populated for this chromosome
            rect(xleft = rep(4, nrow(fdm_windows_merged)), ybottom = fdm_windows_merged$start,
                xright = rep(6, nrow(fdm_windows_merged)), ytop = fdm_windows_merged$end,
                density = -1, col = "#FFF147", border = "#FFF147", lwd = 1.5)
    }
    # plot bars for overlaps
    if(nrow(dxy_fdm_overl) > 0){
        rect(xleft = rep(4, nrow(dxy_fdm_overl)), ybottom = dxy_fdm_overl$start,
            xright = rep(6, nrow(dxy_fdm_overl)), ytop = dxy_fdm_overl$end,
            density = -1, col = "#b5d42a", border = "#b5d42a", lwd = 1.5)
    }
    # plot bars for ELAI regions
    if(nrow(elai_windows_merged) > 0){ # check that table is populated for this chromosome
        rect(xleft = rep(4, nrow(elai_windows_merged)), ybottom = elai_windows_merged$start,
            xright = rep(6, nrow(elai_windows_merged)), ytop = elai_windows_merged$end,
            density = -1, col = "#263CA9", border = "#263CA9", lwd = 1.5)
    }

    # Draw outline of chromosome
    rect(xleft = 4, ybottom = 1, xright = 6, ytop = chr_size, density = 0, 
         border = "black", lwd = 4)
}