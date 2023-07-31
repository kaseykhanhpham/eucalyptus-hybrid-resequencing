# python script to assemble FASTA files of each chloroplast region from Fastplast assemblies
# usage: python make_cp_region_fas.py [sample_name_list] [cp_dir]
# sample_name_list: file with one sample name to include in the tree per line
# cp_dir: directory where sample sequences are located to merge. Must have naming scheme NAME.fas

import os
from sys import argv

sample_name_list = argv[1]
cp_dir = argv[2]

# Read in sample list
sample_list_conn = open(sample_name_list, "r")
sample_list = []
for line in sample_list_conn:
    sample_list.append(line.strip())
sample_list_conn.close()

# Split each cp region into its own file and rename FASTA header to sample ID
for sample in sample_list:
    in_file = "{CP_DIR}/{SAMPLE}.fas".format(CP_DIR = cp_dir, SAMPLE = sample)
    in_conn = open(in_file, "r")
    out_conn = open("dummy_file", "w") # exists so the loop can close the current input file before opening a new one
    for line in in_conn:
        if line.startswith(">"):
            # close previously open file
            out_conn.close()
            # read current chloroplast region
            region = line.strip(">").strip()
            # open relevant region file and write sample ID as header
            out_conn = open("{REGION}_temp.fas".format(REGION = region), "a")
            out_conn.write(">{SAMPLE}\n".format(SAMPLE = sample))
        else:
            # just write line to currently open file if part of sequence body
            out_conn.write(line)
    in_conn.close()
out_conn.close()

os.system("rm dummy_file")