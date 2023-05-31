# python script to calculate Dxy across a VCF file given a size and step for sliding windows and pop structure
# usage: python stats_windows.py [vcf_file] [stat_to_calculate] [window_size] [window_step] [structure_file] [outgroup_file] [chromosomes_list] [out_name]
#
# vcf file: the name/location of the VCF file to perform sliding window stats on
# stat to calculate: must match name of stat in egglib
# window size: size in base pairs of sliding window
# window step: step size in base pairs of sliding window
# structure file: ingroup structure file for egglib in JSON format
#     e.g., {"cluster1":{"eglob_ref":{"samplename":[index], "samplename":[index], etc}, "eglob_MR":
#           {"samplename":[index], "samplename":[index], etc}}}
# outgroup file: outgroup file for egglib in JSON format
# chromosomes list: name of one chromosome from ref genome per line
# out name: name of output file to write

# import libraries
import sys
import egglib
import pandas as pd
import json

# import variables
vcf_name = sys.argv[1]
to_calc = sys.argv[2]
window_size = int(sys.argv[3])
window_step = int(sys.argv[4])
structure_filename = sys.argv[5]
outgroup_filename = sys.argv[6]
chr_list_name = sys.argv[7]
out_name = sys.argv[8]

# import list of chromosomes to iterate through
chr_list = []
chr_list_file = open(chr_list_name, "r")
for line in chr_list_file:
    chr_list.append(line.strip())
chr_list_file.close()

# create 5-col dataframe for storing stat calculations
# 0: chr, 1: start pos, 2: end pos, 3: num variants 4: stat
stats_df = pd.DataFrame(data = None, columns = ["chr", "start", "end", "num_var", to_calc])

# create VCF index
egglib.io.make_vcf_index(vcf_name)

# initiate VCF parser in egglib
vcf_parser = egglib.io.VcfParser(vcf_name)
vcf_parser.load_index()

# create 3-level nested dict to feed to structure object
structure_file = open(structure_filename, "r")
struct_json = json.load(structure_file)
structure_file.close()
outgroup_file = open(outgroup_filename, "r")
outgr_json = json.load(outgroup_file)
outgroup_file.close()

# make ComputeStats object with sliding window stats
pop_struct = egglib.struct_from_dict(struct_json, outgr_json)
computer = egglib.stats.ComputeStats(struct = pop_struct, multi_hits = True)
if to_calc == "FST":
    computer.add_stats("FistWC")
elif to_calc == "Pi":
    computer.add_stats(to_calc)
    computer.add_stats("lseff")
else:
    computer.add_stats(to_calc)

# initiate sliding window object
# iterate through VCF
for chr in chr_list:
    vcf_parser.goto(chr)
    slider = vcf_parser.slider(size = window_size, step = window_step)
    for window in slider:
        # calculate stats from window
        stats_dict = computer.process_sites(window)
        if to_calc == "FST":
            if stats_dict["FistWC"] is None:
                calced_stat = stats_dict["FistWC"]
            else:
                calced_stat = stats_dict["FistWC"][1]
        elif to_calc == "Pi":
            if stats_dict["lseff"] != 0:
                calced_stat = stats_dict[to_calc]/stats_dict["lseff"]
            else:
                calced_stat = stats_dict[to_calc]
        else:
            calced_stat = stats_dict[to_calc]
        # extract position info from window
        chrom = window.chromosome
        start = window.bounds[0]
        end = window.bounds[1]
        num_var = window.num_sites

        # save to dataframe
        stats_df.loc[len(stats_df.index)] = [chrom, start, end, num_var, calced_stat]

# export dataframe as tab-delimited text file
stats_df.to_csv(out_name, sep = "\t", index = False)
