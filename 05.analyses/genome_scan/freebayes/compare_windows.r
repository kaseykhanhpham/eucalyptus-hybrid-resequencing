# R script to compare outlier windows between a set of window tables
# Usage: Rscript compare_windows.r [input_list_1] [input_list_2]
# input_list_1: list of window files to compare against input_list 2
# input_list_2: list of window files to compare against input_list_1
# Comparisons are made pairwise between all possible pairs between input lists 1 and 2
# Dependencies: dplyr library

library(dplyr)

input_list_1 <- read.table(commandArgs(trailingOnly = TRUE)[1], header = FALSE)$V1
input_list_2 <- read.table(commandArgs(trailingOnly = TRUE)[2], header = FALSE)$V1

for(filename1 in input_list_1){
    file1 <- read.table(filename1, header = TRUE)
    for(filename2 in input_list_2){
        if(filename1 != filename2){
            file2 <- read.table(filename2, header = TRUE)
            common_windows <- inner_join(file1, file2)
            outfile_name <- paste("common_", filename1, "_", filename2, sep = "")
            write.table(common_windows, outfile_name, quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
    }
}