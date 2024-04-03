# Rscript to retrieve reciprocal BLAST best matches given two BLAST output files,
# using names of homologous genes from one of the databases 
# Usage: Rscript get_rblast_matches.r -i [blast output 1] -j [blast output 2] 
#                -o [out name] -n [file whose queries are used as names in output]

library(argparser)

# parse arguments
parser <- arg_parser("A script to identify reciprocal best matches from two BLAST output tables.")
parser <- add_argument(parser, "-i", help = "Name of input file 1", type = "character")
parser <- add_argument(parser, "-j", help = "Name of input file 2", type = "character")
parser <- add_argument(parser, "-o", help = "Name of output list", type = "character")
parser <- add_argument(parser, "-n", help = "Number of file (1 or 2) whose query names should be used for output", default = "2", type = "character")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

# import files
infile1 <- read.table(args_list[["i"]], header = FALSE, sep = "\t")
infile2 <- read.table(args_list[["j"]], header = FALSE, sep = "\t")

# get best hit pairs for file 1
infile1_queries <- unique(infile1$V1)
infile1_match <- c()
for(query in infile1_queries){
    subset <- infile1[which(infile1$V1 == query),]
    best_match <- subset[which(subset$V11 == min(subset$V11)),"V2"]
    infile1_match <- c(infile_match, best_match)
}
infile1_best_hits <- data.frame(query = infile1_queries, match = infile1_match)

# get best hit pairs for file 2
infile2_queries <- unique(infile2$V1)
infile2_match <- c()
for(query in infile2_queries){
    subset <- infile2[which(infile2$V1 == query),]
    best_match <- subset[which(subset$V11 == min(subset$V11)),"V2"]
    infile2_match <- c(infile_match, best_match)
}
infile2_best_hits <- data.frame(query = infile2_queries, match = infile2_match)

# Get reciprocal best matches
recip_query <- c()
recip_match <- c()
for(i in c(1:nrow(infile1_best_hits))){
    gene1 <- infile1_best_hits[i, 1]
    gene1_bestmatch <- infile1_best_hits[i, 2]
    # check for presence of pair in file2
    if(gene1_bestmatch %in% infile2_best_hits$V1){
        # check that the bet match has the same genes in file2
        if(infile2_best_hits[which(infile2_bets_hits$V1 == gene1_bestmatch), "V2"] == gene1){
            # add to final list of reciprocal best matches
            if(args_list[["n"]] == "2"){
                recip_query <- c(recip_query, gene1_bestmatch)
                recip_match <- c(recip_match, gene1)
            } else {
                recip_query <- c(recip_query, gene1)
                recip_match <- c(recip_match, gene1_bestmatch)
            }
        }
    }
}

# write final query gene list
write.table(recip_query, args_list[["o"]], quote = FALSE, row.names = FALSE, col.names = FALSE)
