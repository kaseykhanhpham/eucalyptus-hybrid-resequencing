# python script to parse output of BEDTools intersects for gene IDs. Outputs a list of gene IDs.
# Usage: python parse_intersects.py [infile] [outfile]
# infile: BEDTools intersect table with gene annotations retained
# outfile: name of list file to write

import sys

infile_name = sys.argv[1]
outfile_name = sys.argv[2]

gene_id_list = []

infile = open(infile_name, "r")
for line in infile:
    # Parse line by tabs
    cols = line.split("\t")
    # Only extract if column 6 is "gene"
    if cols[5] == "gene":
        # Extract annotation column
        annot = cols[11]
        # Parse annotation
        annot_fields = annot.split(";")
        gene_id = annot_fields[0]
        # Check if field is correct
        if not gene_id.startswith("ID="):
            raise ValueError
        # Store gene ID
        gene_id_list.append("{GENE_ID}\n".format(GENE_ID = gene_id.strip("ID=")))
infile.close()

outfile = open(outfile_name, "w")
outfile.writelines(gene_id_list)
outfile.close()