# R script to interpolate genetic distances between variant positions
# in a counts file from vcftools --count provided a genetic map
# Usage: Rscript get_gen_dists.r [counts_file] [genetic_map] [out_name]
# counts_file: the output of vcftools --counts
# genetic_map: a genetic map with the columns chr, variant ID, genetic distance, physical distance in that order.
# out_name: name of the output file to be written with distances corresponding to each row of the provided counts table

# turn off scientific notation
options(scipen=999)

# import CLI arguments
counts_filename <- commandArgs(trailingOnly = TRUE)[1]
map_filename <- commandArgs(trailingOnly = TRUE)[2]
out_name <- commandArgs(trailingOnly = TRUE)[3]

counts_tab <- read.table(counts_filename, skip = 1, header = FALSE, col.names = c("Chr", "Pos", "nalleles", "totalalleles", "A_count", "a_count"), colClasses = c("character", "numeric", "numeric", "numeric", "character", "character"))
map_tab <- read.table(map_filename, header = FALSE, col.names = c("Chr", "ID", "cM", "bp"), colClasses = c("character", "character", "numeric", "numeric"))

# define function to interpolate genetic distance if passed a chromosome, position, and genetic map.
# Debug test cases:
# Chr03 53000000, Expected: 8640.1
# Chr06 50600000, Expected: 13799.26
# Chr05 969243, Expected: 0.00
# Chr05 500000, Expected: -0.4747846
# Chr09 39043096, Expected: 8628.917
get_dist <- function(chr, pos, genmap){
    write(paste("doing", chr, pos, sep = " "), stdout()) # debug

    all_chr_pos <- genmap[which(genmap$Chr == chr), "bp"]
    smaller <- which(all_chr_pos <= pos)
    bigger <- which(all_chr_pos >= pos)

    # account for edgecases (literally)
    if(length(smaller) == 0){
        # calculate cM per bp for first chunk of map
        first_pos <- genmap[which(genmap$Chr == chr), "bp"][1]
        first_dist <- genmap[which(genmap$Chr == chr), "cM"][1]
        second_pos <- genmap[which(genmap$Chr == chr), "bp"][2]
        second_dist <- genmap[which(genmap$Chr == chr), "cM"][2]
        # if adjacent map markers have same genetic dist (i.e., rate of cM per bp will be 0),
        # walk forwards until a marker with a different distance is found
        step_counter = 3
        while(second_dist == first_dist){
            second_pos <- genmap[which(genmap$Chr == chr), "bp"][step_counter]
            second_dist <- genmap[which(genmap$Chr == chr), "cM"][step_counter]
            step_counter <- step_counter + 1
        }
        cm_per_bp <- (second_dist - first_dist)/(second_pos - first_pos)
        # multiply by how far back the passed position is
        dist <- first_dist - cm_per_bp*(first_pos - pos)

    } else if(length(bigger) == 0){
        # calculate cM per bp for the last chunk of the map
        last_i <- length(which(genmap$Chr == chr))
        ult_pos <- genmap[which(genmap$Chr == chr), "bp"][last_i]
        ult_dist <- genmap[which(genmap$Chr == chr), "cM"][last_i]
        penult_pos <- genmap[which(genmap$Chr == chr), "bp"][last_i - 1]
        penult_dist <- genmap[which(genmap$Chr == chr), "cM"][last_i - 1]
        # if adjacent map markers have same genetic dist (i.e., rate of cM per bp will be 0),
        # walk backwards until a marker with a different distance is found
        step_counter = last_i - 2
        while(penult_dist == ult_dist){
            penult_pos <- genmap[which(genmap$Chr == chr), "bp"][step_counter]
            penult_dist <- genmap[which(genmap$Chr == chr), "cM"][step_counter]
            step_counter <- step_counter - 1
        }
        cm_per_bp <- (ult_dist - penult_dist)/(ult_pos - penult_pos)
        # multiply by how far ahead the passed position is
        dist <- ult_dist + cm_per_bp*(pos - ult_pos)

    } else {
        # variant is not before the first recorded distance or after the last,
        # retrieve its bracketing positions
        before_i <- max(smaller)
        after_i <- min(bigger)

        # check for if variant's location is actually precisely in genetic map
        if(before_i == after_i){
            dist <- genmap[which(genmap$Chr == chr), "cM"][before_i]
        } else {
            before_dist <- genmap[which(genmap$Chr == chr), "cM"][before_i]
            after_dist <- genmap[which(genmap$Chr == chr), "cM"][after_i]

            before_pos <- genmap[which(genmap$Chr == chr), "bp"][before_i]
            after_pos <- genmap[which(genmap$Chr == chr), "bp"][after_i]

            fraction_between <- ((pos - before_pos)/(after_pos - before_pos))
            dist <- before_dist + fraction_between * (after_dist - before_dist)
        }
    }
    morgans <- dist*100
    return(morgans)
}

absolute_dists <- apply(counts_tab, 1, function(curr_row) get_dist(curr_row[1], as.numeric(curr_row[2]), map_tab))
rel_dists <- sapply(c(2:length(absolute_dists)), function(x) absolute_dists[x] - absolute_dists[x - 1])
rel_dists <- c(0, rel_dists)

out_tab <- data.frame(cbind(counts_tab$Chr, counts_tab$Pos, rel_dists))

write.table(out_tab, out_name, quote = FALSE, row.names = FALSE, col.names = FALSE)