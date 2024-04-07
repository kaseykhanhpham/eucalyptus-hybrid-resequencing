# Rscript to retrieve reciprocal BLAST best matches given two BLAST output files,
# using names of homologous genes from one of the databases 
# Usage: Rscript get_rblast_matches.r -i [blast output 1] -j [blast output 2] 
#                -l [out_list_name] -t [out table_name] -n [file whose queries are used as names in output]

library(argparser)

# parse arguments
parser <- arg_parser("A script to identify reciprocal best matches from two BLAST output tables.")
parser <- add_argument(parser, "-i", help = "Name of input file 1", type = "character")
parser <- add_argument(parser, "-j", help = "Name of input file 2", type = "character")
parser <- add_argument(parser, "-l", help = "Name of output list", type = "character")
parser <- add_argument(parser, "-t", help = "Name of output table", type = "character")
parser <- add_argument(parser, "-n", help = "Number of file (1 or 2) whose query names should be used for output", default = "2", type = "character")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

# import files
infile1 <- read.table(args_list[["i"]], header = FALSE, sep = "\t")
infile2 <- read.table(args_list[["j"]], header = FALSE, sep = "\t")

# get best hit pairs for file 1
infile1_queries <- unique(infile1$V1)
infile1_matches <- list()
for(query in infile1_queries){
    subset <- infile1[which(infile1$V1 == query),]
    best_match <- subset[which(subset$V11 == min(subset$V11)),"V2"]
    infile1_matches[[query]] <- best_match
}

# get best hit pairs for file 2
infile2_queries <- unique(infile2$V1)
infile2_matches <- list()
for(query in infile2_queries){
    subset <- infile2[which(infile2$V1 == query),]
    best_match <- subset[which(subset$V11 == min(subset$V11)),"V2"]
    infile2_matches[[query]] <- best_match
}

# Get reciprocal best matches
recip_query <- c()
recip_match <- c()
for(gene1 in names(infile1_matches)){
    gene1_bestmatch <- infile1_matches[[gene1]]
    # check for presence of pair in file2
    # single hit for gene1
    if(length(gene1_bestmatch) == 1){
        # check if reverse match is present in file2
        gene2 <- gene1_bestmatch
        recip_bool <- gene1 %in% infile2_matches[[gene1_bestmatch]]
        # note that this allows there to be mult best matches in file2 as long as gene1 is one of them.
    # mult hits for gene1
    } else {
        # check for how many of the matches gene1 is in the reverse file2 hits
        recip_hits <- c()
        for(match in gene1_bestmatch){
            # print(infile2_matches[[match]]) # DEBUG
            hit_bool <- gene1 %in% infile2_matches[[match]]
            # print(hit_bool) # DEBUG
            recip_hits <- c(recip_hits, hit_bool)
        }
        # if no match...
        if(length(which(recip_hits)) == 0){
            gene2 <- NA
            recip_bool <- FALSE
        # if only one match...
        } else if(length(which(recip_hits)) == 1){
            gene2 <- gene1_bestmatch[which(recip_hits)]
            recip_bool <- TRUE
        # if multiple matches...
        } else {
            # pick one of the matches at random
            hit_to_use <- sample(c(1:length(which(recip_hits))), 1, replace = FALSE)
            gene2 <- gene1_bestmatch[which(recip_hits)[hit_to_use]]
            recip_bool <- TRUE
        }
    }
    # record final reciprocal best match from selection process
    if(recip_bool){
        if(args_list[["n"]] == "2"){
            recip_query <- c(recip_query, gene2)
            recip_match <- c(recip_match, gene1)
        } else {
            recip_query <- c(recip_query, gene1)
            recip_match <- c(recip_match, gene2)
        }
    }
}

# write final query gene list
recip_best_hits <- data.frame(query = recip_query, match = recip_match)
write.table(recip_best_hits, args_list[["t"]], quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
recip_query_prot <- paste(recip_query, ".p", sep = "")
write.table(recip_query_prot, args_list[["l"]], quote = FALSE, row.names = FALSE, col.names = FALSE)
