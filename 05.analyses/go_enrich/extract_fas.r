# R script to retrieve FASTA file given a source FASTA file and list of headers to extract

library(argparser)
library(seqinr)

# parse arguments
parser <- arg_parser("A script to extract a subset of FASTA sequences given a list of headers.")
parser <- add_argument(parser, "-i", help = "Name of FASTA input", type = "character")
parser <- add_argument(parser, "-l", help = "List of headers to extract", type = "character")
parser <- add_argument(parser, "-o", help = "Name of output FASTA to write", type = "character")
parser <- add_argument(parser, "-t", help = "Sequence type (DNA or AA)", type = "character", default = "AA")

args_list <- parse_args(parser, argv = commandArgs(trailingOnly = TRUE))

infas <- read.fasta(file = args_list[["i"]], seqtype = args_list[["t"]])
seqlist <- read.table(args_list[["l"]], header = FALSE)$V1

# get index of sequences in list
out_inds <- match(seqlist, names(infas))
# subset FASTA file
outfas <- infas[seqlist]
# output subsetted FASTA file
write.fasta(outfas, names(outfas), args_list[["o"]])
