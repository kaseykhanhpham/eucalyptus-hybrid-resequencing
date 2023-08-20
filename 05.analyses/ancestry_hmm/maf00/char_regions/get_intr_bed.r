# R script to extract introgressed regions from the output posterior file of Ancestry_HMM
# Usage: get_intr_bed.r [posterior_file] [seed_posterior] [minimum_posterior] [include_dist] [warning_posterior] [output_name]
#        posterior_file: the output file of AncestryHMM for a single sample
#                        (assumes single pulse of introgression in a two-donor system)
#        seed_posterior: at least one variant within an identified region must have at least this posterior value
#        minimum_posterior: variants more than [include_dist] apart with posteriors below this value signal the end of the identified region
#        include_dist: in bp, see minimum_posterior
#        warning_posterior: a warning will be printed upon inclusion of variants under this value
#        output_name: name of the BED file to be printed with final regions

options(scipen = 999)

filename <- commandArgs(trailingOnly = TRUE)[1]
seed_post <- as.numeric(commandArgs(trailingOnly = TRUE)[2])
min_post <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
include_dist <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
warn_post <- as.numeric(commandArgs(trailingOnly = TRUE)[5])
output_name <- commandArgs(trailingOnly = TRUE)[6]

infile_colnames <- c("chrom", "pos", "hom_glob", "het", "hom_cord")
post_tab <- read.table(filename, skip = 1, sep = "\t", col.names = infile_colnames, colClasses = c("character", "integer", "numeric", "numeric", "numeric"))
chrs <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

chr_regions <- list()

for(chr in chrs){
    region_seeds <- which(post_tab[,5] >= seed_post)
    min_regions <- which(post_tab[,5] >= min_post)
    # subset to just the current chromosome
    rows_in_chr <- which(post_tab[,1] == chr)
    region_seeds <- intersect(region_seeds, rows_in_chr)
    min_regions <- intersect(min_regions, rows_in_chr)

    regions_present <- (length(min_regions) > 0) & (length(region_seeds) > 0)
    if(regions_present){
        # merge region fragments separated by less than 5 variants
        gaps_to_fill <- unlist(sapply(c(2:length(min_regions)), function(i) {
            curr <- min_regions[i]
            prev <- min_regions[i - 1]

            curr_pos <- post_tab[curr,2]
            prev_pos <- post_tab[prev,2]

            need_to_fill <- (curr - prev > 1) & (curr_pos - prev_pos < include_dist)
            if(need_to_fill){
                return(c((prev + 1):(curr - 1)))
            } else {
                return(NA)
            }
        }))
        gaps_to_fill <- gaps_to_fill[which(!is.na(gaps_to_fill))]
        min_regions <- sort(c(min_regions, gaps_to_fill))

        # Warn user if included a variant below 85% posterior probability
        min_region_posts <- post_tab[min_regions, 5]
        low_posts_pos <- post_tab[min_regions, 2][which(min_region_posts < warn_post)]
        if(any(min_region_posts < warn_post)){
            write(paste("WARNING: Included variants at", chr, "position", low_posts_pos, "below", warn_post, "posterior probability in region."), stdout())
        }

        # get starts and ends of filled-out minimum posterior regions
        # first, get indices of start and end rows
        start_indices <- unlist(sapply(c(2:length(min_regions)), function(i) {
            curr <- min_regions[i]
            prev <- min_regions[i - 1]
            if((curr - prev) > 1){
                return(curr)
            } else {
                return(NA)
            }
        }))
        start_indices <- start_indices[which(!is.na(start_indices))]
        start_indices <- c(min_regions[1], start_indices)

        end_indices <- unlist(sapply(c(2:length(min_regions)), function(i) {
            curr <- min_regions[i]
            prev <- min_regions[i - 1]
            if((curr - prev) > 1){
                return(prev)
            } else {
                return(NA)
            }
        }))
        end_indices <- end_indices[which(!is.na(end_indices))]
        end_indices <- c(end_indices, min_regions[length(min_regions)])

        # then, extract actual bp positions from indices
        if(length(start_indices) != length(end_indices)){
            write(paste("WARNING: number of region starts does not equal number of region ends for", chr), stdout())
        }
        reg_starts <- post_tab[start_indices, 2]
        reg_ends <- post_tab[end_indices, 2]

        # filter minimum regions to those that contain seeds
        reg_with_seeds_indices <- c()
        for(seed in region_seeds){
            bigger_than <- which(reg_starts <= post_tab[seed,2])
            smaller_than <- which(reg_ends >= post_tab[seed,2])
            withindex <- intersect(bigger_than, smaller_than) # should only be one
            reg_with_seeds_indices <- c(reg_with_seeds_indices, withindex)
        }
        # get rid of any duplicates from multiple seeds in the same region
        reg_with_seeds_indices <- unique(reg_with_seeds_indices)

        # Make df of final regions for the chromosome
        reg_starts_final <- reg_starts[reg_with_seeds_indices]
        reg_ends_final <- reg_ends[reg_with_seeds_indices]
        chr_col <- rep(chr, length(reg_starts_final))
        regions_df <- as.data.frame(cbind(chr = chr_col,
                                        start = reg_starts_final,
                                        end = reg_ends_final))
    } else {
        # if no regions above .90/.95 posterior probability, store empty dataframe
        regions_df <- data.frame(chr = character(), start = numeric(), end = numeric())
    }
    # Store in global list under the chromosome name
    chr_regions[[chr]] <- regions_df
}

# concatenate dfs for each chromosome
bed_df <- chr_regions[[chrs[1]]]
for(i in c(2:length(chrs))){
    bed_df <- rbind(bed_df, chr_regions[[chrs[i]]])
}

# print to BED file
# (I actually don't know whether AncestryHMM lists positions with 0- or
#  1-indexing, so this may be one bp off according to BED format. Oh well.)
write.table(bed_df, output_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
