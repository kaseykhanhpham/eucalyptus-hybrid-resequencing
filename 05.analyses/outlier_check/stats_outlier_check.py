# python script to calculate dXY pairwise for outlier WF03 against other populations on average
# One-off thing; files hard-coded into script.

# import libraries
import sys
import egglib
import pandas as pd
import json
import copy
 
# files to read in
vcf_name = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
structure_filename = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/glob_structure.json"
outgroup_filename = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/glob_outgroup.json"
chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "ChrUn"]
out_root = "stats_outlier_check"

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
outgroup_file = open(outgroup_filename, "r") # outgroup is E. cordata + SRR10339635
outgr_json = json.load(outgroup_file)
outgroup_file.close()

# iterate through all reference samples of E. globulus, excluding each and then calculating average pi over the genome
for key in struct_json["cluster1"]["eglob_ref"]:
    # make deep copies of base population structure dicts
    struct_json_copy = copy.deepcopy(struct_json)
    outgr_json_copy = copy.deepcopy(outgr_json)

    # edit deep copies to exclude E. glob reference sample
    rem_val = struct_json_copy["cluster1"]["eglob_ref"][key]
    del struct_json_copy["cluster1"]["eglob_ref"][key] # remove from ingroup structure
    outgr_json_copy[key] = rem_val # add to outgroup dict

    # initialize egglib computer for genome stat calculations and dataframe to store calcs
    pop_struct = egglib.struct_from_dict(struct_json_copy, outgr_json_copy)
    computer = egglib.stats.ComputeStats(struct = pop_struct, multi_hits = True)
    computer.add_stats("Pi", "Dxy")
    stats_df = pd.DataFrame(data = None, columns = ["chr", "start", "end", "Pi", "Dxy"])
    
    # initialize sliding window object and iterate through each chromosome
    for chr in chr_list:
        vcf_parser.goto(chr)
        slider = vcf_parser.slider(size = 5000, step = 2500)
        for window in slider:
            # calculate stats from window
            stats_dict = computer.process_sites(window)
            pi = stats_dict["Pi"]
            dxy = stats_dict["Dxy"]
            # extract position info from window
            chrom = window.chromosome
            start = window.bounds[0]
            end = window.bounds[1]
            # save to dataframe
            stats_df.loc[len(stats_df.index)] = [chrom, start, end, pi, dxy]
    
    # Write to output file
    out_name = "{OUT_ROOT}_no_{KEY}.tab".format(OUT_ROOT = out_root, KEY = key)
    stats_df.to_csv(out_name, sep = "\t", index = False)

# And then plot resulting distribution of Pi and Dxy to see whether excluding WF03 makes a difference!