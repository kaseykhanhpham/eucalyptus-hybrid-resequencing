# R script to reformat genome annotation in GFF format to BED12 format with defined exons in blocks
# Usage: Rscript get_annot_bed.r [infile] [outfile]
#     infile: input GFF annotation file
#     outfile: name of output BED12 file to write
# BED12: chr start end name score strand thickstart thickend itemRbg blockcount blocksizes blockstarts
# GFF: seqname source feature start end score strand frame attr

infile_name <- commandArgs(trailingOnly = TRUE)[1]
outfile_name <- commandArgs(trailingOnly = TRUE)[2]

infile <- read.table(infile_name, header = FALSE, sep = "\t")

just_genes <- infile[which(infile$V3 == "gene"),]
just_cds <- infile[which(infile$V3 == "CDS"),]

# get unique gene names
gene_names <- sapply(just_genes$V9, function(x) unlist(strsplit(x, ";", fixed = TRUE))[1], USE.NAMES = FALSE)
gene_names <- gsub("ID=", "", gene_names)

# get overall gene of each CDS
gene_membership <- sapply(just_cds$V9, function(x) unlist(strsplit(x, ";", fixed = TRUE))[1], USE.NAMES = FALSE)
gene_membership <- gsub("ID=", "", gene_membership)
gene_membership <- gsub("-RA:cds", "", gene_membership)

# remove any genes in the gene-only table and list of gene names that don't
# have at least one associated CDS annotation
has_cds <- which(gene_names %in% gene_membership)
just_genes <- just_genes[has_cds,]
gene_names <- gene_names[has_cds]

# List of lists of indices for rows that belong to each gene
gene_pos <- sapply(gene_names, function(gene) which(gene_membership == gene))

# Get block count, block sizes, and block starts
blockcount <- sapply(gene_pos, function(pos_vec) length(pos_vec))
block_sizes <- sapply(gene_pos, function(pos_vec) 
    sapply(pos_vec, function(pos) 
        just_cds[pos, "V5"] - just_cds[pos, "V4"]+ 1))
block_start <- sapply(gene_pos, function(pos_vec) 
    sapply(pos_vec, function(pos) 
        just_cds[pos, "V4"] - 1))

# construct final BED12 dataframe
chr <- just_genes$V1
start <- just_genes$V4 - 1
end <- just_genes$V5
name <- gene_names
score <- rep(0, nrow(just_genes))
strand <- just_genes$V7
thickstart <- rep(0, nrow(just_genes))
thickend <- rep(0, nrow(just_genes))
itemRbg <- rep(0, nrow(just_genes))
# blockcounts is fine the way it is
blocksizes <- sapply(block_sizes, function(x) paste(x, collapse=","))
blockstarts <- sapply(block_start, function(x) paste(x, collapse=","))

final_bed <- cbind(chr, start, end, name, score, strand, 
thickstart, thickend, itemRbg, 
blockcount, blocksizes, blockstarts)

write.table(final_bed, outfile_name, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")