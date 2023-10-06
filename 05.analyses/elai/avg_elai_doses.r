# R script to plot the results of ELAI for an individual
# Usage: R avg_elai_doses.r [files_to_average] [output_name]
#     files_to_average: a list of ps21 files from ELAI output over which to average,
#                       one per line, full and relative addreses from CWD OK
#     output_name: the name of the output file to write

options(scipen = 999)

infile_list_name <- commandArgs(trailingOnly = TRUE)[1]
outname <- commandArgs(trailingOnly = TRUE)[2]

# import input files from cmd line
infile_list <- read.table(infile_list_name, header = FALSE, as.is = TRUE)$V1

# read in all dosage files into a list
p21_list <- list()
for(i in c(1:length(infile_list))){
    p21_list[[i]] <- read.table(infile_list[i], header = FALSE)
}

# get dimensions of matrix to print (should be the same as for all inputs)
num_rows <- nrow(p21_list[[1]])
num_cols <- ncol(p21_list[[1]])

# initialize empty output matrix
avg_mat <- matrix(, nrow = num_rows, ncol = num_cols)

# iterate through positions in matrix, averaging across input list
for(i in c(1:num_rows)){
    for(j in c(1:num_cols)){
        cell_vals <- unlist(lapply(p21_list, function(x) x[i,j]))
        cell_avg <- mean(cell_vals)
        avg_mat[i,j] <- round(cell_avg, 3)
    }
}

write.table(avg_mat, outname, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
write(paste("Done writing to", outname), stdout())
