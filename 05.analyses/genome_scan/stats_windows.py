# python script to calculate Dxy across a VCF file given a size and step for sliding windows and pop structure
# usage: python stats_windows.py [vcf_file] [window_size] [window_step] [structure_file] [outgroup_file] [chromosomes_list] [out_name]

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
outgroup_filename = sys.argv[5]
chr_list_name = sys.argv[6]
out_name = sys.argv[7]

# import list of chromosomes to iterate through
chr_list = []
chr_list_file = open(chr_list_name, "r")
for line in chr_list_file:
    chr_list.append(line.strip())
chr_list_file.close()

# create 4-col dataframe for storing Dxy calculations
# 0: chr, 1: start pos, 2: end pos, 3: avg Dxy
stats_df = pd.DataFrame(data = None, columns = ["chr", "start", "end", "Pi", "Deta", "Dxy"])

# create VCF index
egglib.io.make_vcf_index(vcf_name)

# initiate VCF parser in egglib
vcf_parser = egglib.io.VcfParser(vcf_name)
vcf_parser.load_index()

# create 3-level nested dict to feed to structure object
# struct_dict = {clust1:{glob_ref:{samplename:[], samplename:[], etc}, glob_MR:{samplename:[], samplename:[], etc}}}
structure_file = open(structure_filename, "r")
struct_json = json.load(structure_file)
structure_file.close()
outgroup_file = open(outgroup_filename, "r")
outgr_json = json.load(outgroup_file)
outgroup_file.close()

# make ComputeStats object with sliding window stats
pop_struct = egglib.struct_from_dict(struct_json, outgr_json)
computer = egglib.stats.ComputeStats(struct = pop_struct, multi_hits = True)
computer.add_stats("Pi", "Deta", "Dxy")

# initialize storage table
stats_df = pd.DataFrame(data = None, columns = ["chr", "start", "end", "Pi", "Deta", "Dxy"])

# initiate sliding window object
# iterate through VCF
for chr in chr_list:
    vcf_parser.goto(chr)
    slider = vcf_parser.slider(size = window_size, step = window_step)
    for window in slider:
        # calculate stats from window
        stats_dict = computer.process_sites(window)
        pi = stats_dict["Pi"]
        deta = stats_dict["Deta"]
        dxy = stats_dict["Dxy"]
        # extract position info from window
        chrom = window.chromosome
        start = window.bounds[0]
        end = window.bounds[1]
        # save to dataframe
        stats_df.loc[len(stats_df.index)] = [chrom, start, end, pi, deta, dxy]

# export dataframe as tab-delimited text file
stats_df.to_csv(out_name, sep = "\t", index = False)