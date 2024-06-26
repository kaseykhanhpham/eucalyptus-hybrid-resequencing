# Set of functions to use for quickly generating heatmaps to visualize
# the output of AncestryHMM

get_windows <- function(chr_name, start, end, num_splits){
    total_size <- end - start + 1
    split_size <- round(total_size/num_splits)
    starts <- c(start)
    ends <- c()
    for(i in c(1:(num_splits - 1))){
        starts <- c(starts, (start + ((i) * split_size)))
        ends <- c(ends, (start + ((i) * split_size) - 1))
    }
    ends <- c(ends, end)
    chr <- rep(chr_name, num_splits)
    out_tab <- data.frame(chrom = chr, start = starts, end = ends)
    return(out_tab)
}

populate_mat_val <- function(sample_vec, window_tab, post_file_loc){
    mat_list <- list()
    index = 1
    window_order <- c()
    for(sample in sample_vec){
        infile_name <- paste(post_file_loc, "/", sample, ".posterior", sep = "")
        infile_colnames <- c("chrom", "pos", "hom_glob", "het", "hom_cord")
        infile <- read.table(infile_name, skip = 1, sep = "\t", col.names = infile_colnames, colClasses = c("character", "integer", "numeric", "numeric", "numeric"))
        # construct row for each window -- 1 for a cell if the ind has any homozygous cord posterior > cutoff in the window
        # 0.5 for a cell if the ind has any heterozygous cord posterior > cutoff in the window
        for(i in c(1:nrow(window_tab))) {
            wchr <- window_tab[i, "chrom"]
            wsta <- window_tab[i, "start"]
            wend <- window_tab[i, "end"]
            infile_mask <- which(infile$chrom == wchr & infile$pos >= wsta & infile$pos <= wend)
            infile_subset <- infile[infile_mask,]
            hom_cell <- max(infile_subset$hom_cord)
            het_cell <- max(infile_subset$het)
            
            windex <- paste(wchr, wsta, wend, sep = "_")
            window_order <- c(window_order, windex)
            mat_ind <- paste("row", index, sep = "_")

            mat_list[[mat_ind]]["window"] <- windex
            mat_list[[mat_ind]]["sample"] <- sample
            mat_list[[mat_ind]]["hom_mat"] <- hom_cell
            mat_list[[mat_ind]]["het_mat"] <- het_cell
            
            index <- index + 1
        }
    }
    mat_df <- as.data.frame(t(as.data.frame(mat_list)))
    mat_df$hom_mat <- as.numeric(mat_df$hom_mat)
    mat_df$het_mat <- as.numeric(mat_df$het_mat)
    # convert window indexes to factor so they will appear in the same order as the input table
    window_order <- unique(window_order)
    mat_df$window <- factor(mat_df$window, levels = window_order)
    return(mat_df)
}

plot_heatmap <- function(heat_mat, names_from, names_to, chr_name, title_epithet){
    # heat_mat: the output of populate_mat()
    # names_from: a vector of sample names to translate from in the heatmap matrix
    # names_to: a vector of sample names to translate to in final plots

    library(ggplot2)

    label_order <- match(heat_mat$sample, names_from)
    heat_mat$acc <- names_to[label_order]
    plot_xlabs <- unique(gsub(paste(chr_name, "_", sep = ""), "", heat_mat$window))

    hom_plot <- ggplot(heat_mat, aes(window, acc, fill = hom_mat)) + 
                geom_tile(color = "#004c5f") + 
                scale_fill_steps(high = "#ffd966", low = "#0087a8") + 
                guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
                ggtitle(paste("Largest posterior of homozygous E. cordata ancestry", title_epithet)) + 
                xlab(paste(chr_name, "Window")) + 
                ylab("E. globulus Sample") + 
                scale_x_discrete(labels = plot_xlabs) +
                theme(plot.title=element_text(size=16),
                  axis.title = element_text(size = 14),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 14))

    het_plot <- ggplot(heat_mat, aes(window, acc, fill = het_mat)) + 
                geom_tile(color = "#004c5f") + 
                scale_fill_steps(high = "#ffd966", low = "#0087a8") + 
                guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
                ggtitle(paste("Largest posterior of heterozygous E. cordata ancestry", title_epithet)) + 
                xlab(paste(chr_name, "Window")) + 
                ylab("E. globulus Sample") + 
                scale_x_discrete(labels = plot_xlabs) +
                theme(plot.title=element_text(size=16),
                  axis.title = element_text(size = 14),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 14))
    return(list(hom_plot = hom_plot, het_plot = het_plot))
}