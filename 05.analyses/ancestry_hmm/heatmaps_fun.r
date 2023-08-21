# Set of functions to use for quickly generating heatmaps to visualize
# the output of AncestryHMM

populate_mat <- function(sample_vec, window_tab, post_cutoff){
    mat_list <- list()
    index = 1
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
            if(any(infile_subset$hom_cord >= post_cutoff)){
                hom_cell <- 1
                het_cell <- 1
            } else if(any(infile_subset$het >= post_cutoff)){
                hom_cell <- 0
                het_cell <- 0.5
            } else {
                hom_cell <- 0
                het_cell <- 0
            }
            windex <- paste(wchr, wsta, wend, sep = "_")
            mat_ind <- paste("row", index, sep = "_")

            mat_list[[mat_ind]]["window"] <- windex
            mat_list[[mat_ind]]["sample"] <- sample
            mat_list[[mat_ind]]["hom_mat"] <- hom_cell
            mat_list[[mat_ind]]["het_mat"] <- het_cell
            
            index <- index + 1
        }
    }
    mat_df <- as.data.frame(t(as.data.frame(mat_list)))
    mat_df$hom_mat <- as.integer(mat_df$hom_mat)
    mat_df$het_mat <- as.numeric(mat_df$het_mat)
    return(mat_df)
}

plot_heatmap <- function(heat_mat, names_from, names_to, chr_name, post_cutoff){
    # heat_mat: the output of populate_mat()
    # names_from: a vector of sample names to translate from in the heatmap matrix
    # names_to: a vector of sample names to translate to in final plots

    library(ggplot2)

    label_order <- match(heat_mat$sample, names_from)
    heat_mat$acc <- names_to[label_order]
    plot_xlabs <- unique(gsub(chr_name, "", heat_mat$window))

    hom_plot <- ggplot(heat_mat, aes(window, acc, fill = hom_mat)) + 
                geom_tile(color = "#004c5f") + 
                scale_fill_steps(high = "#ffd966", low = "#0087a8") + 
                guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
                ggtitle(paste("At least one locus with >", post_cutoff, "homozygous E. cordata ancestry")) + 
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
                ggtitle(paste("At least one locus with >", post_cutoff, "E. cordata ancestry")) + 
                xlab(paste(chr_name, "Window")) + 
                ylab("E. globulus Sample") + 
                scale_x_discrete(labels = plot_xlabs) +
                theme(plot.title=element_text(size=16),
                  axis.title = element_text(size = 14),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 14))
    return(list(hom_plot = hom_plot, het_plot = het_plot))
}

