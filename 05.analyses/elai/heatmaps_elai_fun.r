# Set of functions to use for quickly generating heatmaps to visualize
# the output of Easy Local Ancestry Inference (ELAI)

populate_elai_mat_val <- function(infile_loc, sample_vec, window_tab){
    mat_list <- list()
    index = 1
    for(i in c(1:nrow(window_tab))) {
        wchr <- window_tab[i, "chrom"]
        wsta <- window_tab[i, "start"]
        wend <- window_tab[i, "end"]
        # import the corresponding chromosome files
        dosage_name <- paste(infile_loc, "/", wchr, "_avg.ps21.txt", sep = "")
        pos_name <- paste(infile_loc, "/", wchr, "_avg.snpinfo.txt", sep = "")
        dosage_in <- scan(dosage_name)
        pos_in <- read.table(pos_name, header = TRUE, sep = "\t")
        dim(dosage_in) <- c(2,nrow(pos_in),20) # 2 ancestral populations, num snps, 20 inds

        ## subset to just what is in the current window
        # get just chromosome for each position in the snpinfo file
        pos_in_chr <- unlist(lapply(strsplit(pos_in$rs, ":"), function(x) x[1]))
        # subset dosage file to just E. cordata ancestry in current window
        infile_mask <- which(pos_in$pos >= wsta & pos_in$pos <= wend)
        infile_subset <- dosage_in[2,infile_mask,]

        # iterate through each sample and add row in matrix list for each
        for(i in c(1:length(sample_vec))){
            # get maximum dosage in window for a sample
            w_max_dose <- max(infile_subset[,i])
            # make indices for output matrix
            windex <- paste(wchr, wsta, wend, sep = "_")
            mat_ind <- paste("row", index, sep = "_")
            # populate output matrix
            mat_list[[mat_ind]]["window"] <- windex
            mat_list[[mat_ind]]["sample"] <- sample_vec[i]
            mat_list[[mat_ind]]["dosage"] <- w_max_dose
        
            index <- index + 1
        }
    }
    mat_df <- as.data.frame(t(as.data.frame(mat_list)))
    mat_df$dosage <- as.numeric(mat_df$dosage)
    return(mat_df)
}

plot_elai_heatmap <- function(heat_mat, names_from, names_to, chr_name){
    # heat_mat: the output of populate_mat()
    # names_from: a vector of sample names to translate from in the heatmap matrix
    # names_to: a vector of sample names to translate to in final plots

    library(ggplot2)

    label_order <- match(heat_mat$sample, names_from)
    heat_mat$acc <- names_to[label_order]
    plot_xlabs <- unique(gsub(chr_name, "", heat_mat$window))

    dose_plot <- ggplot(heat_mat, aes(window, acc, fill = dosage)) + 
                geom_tile(color = "#000000") + 
                scale_fill_steps(high = "#616161", low = "#ffffff") + 
                guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
                ggtitle(paste(chr_name, "maximum dosage of E. cordata ancestry per window")) + 
                xlab(paste(chr_name, "Window")) + 
                ylab("E. globulus Sample") + 
                scale_x_discrete(labels = plot_xlabs) +
                theme(plot.title=element_text(size=16),
                      axis.title = element_text(size = 14),
                      axis.text.x = element_text(size = 9),
                      axis.text.y = element_text(size = 14))
}

