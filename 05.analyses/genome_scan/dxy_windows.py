# python script to calculate Dxy across a VCF file given a size and step for sliding windows and pop structure
# usage: python dxy_windows.py [vcf_file] [window_size] [window_step] [structure_file]

# import libraries
import sys
import egglib
import pandas as pd
import json

# import variables
vcf_name = sys.argv[1]
window_size = int(sys.argv[2])
window_step = int(sys.argv[3])
structure_filename = sys.argv[4]

# create 4-col dataframe for storing Dxy calculations
# 0: chr, 1: start pos, 2: end pos, 3: avg Dxy
dxy_df = pd.DataFrame(data = None, columns = ["chr", "start", "end", "Dxy"])

# initiate VCF parser in egglib
vcf_parser = egglib.io.VcfParser(vcf_name)

# create 3-level nested dict to feed to structure object
# struct_dict = {clust1:{glob_ref:{samplename:[], samplename:[], etc}, glob_MR:{samplename:[], samplename:[], etc}}}
structure_file = open(structure_filename, "r")
struct_json = json.load(structure_file)
structure_file.close()
pop_struct = egglib.struct_from_dict(struct_json, None)

# create ComputeStats instance for Dxy
computer = egglib.stats.ComputeStats(struct = pop_struct, multi_hits = True)
computer.add_stats("Pi", "Deta", "Dxy")

# initiate sliding window object
# window: 100000, step: 20000
# iterate through VCF
slider = vcf_parser.slider(size = 100000, step = 20000)
for window in slider:
    # calculate Dxy from window
    dxy = computer.process_sites(window).values()[0] # EXTRACT FROM DICT... CHECK THIS
    # extract position info from window
    chrom = window.chromosome
    start = window.bounds[0]
    end = window.bounds[1]
    # save to dataframe
    dxy_df.loc[len(dxy_df.index)] = [chrom, start, end, dxy]

# export dataframe as tab-delimited text file
dxy_df.to_csv("{INFILE}_Dxy.tab".format(INFILE = vcf_name), sep = "\t")
