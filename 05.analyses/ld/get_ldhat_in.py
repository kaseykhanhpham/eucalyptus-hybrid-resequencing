# python script to format VCF file for input into LDhat using genotypes (unphased data)
# assumes VCF used is sorted, already filtered to biallelic SNPs only and contains only samples of interest
# meant to be run separately for each chromosome in a genome
# outputs a separate file for each chromosome listed in the VCF file
# Usage: python get_ldhat_in.py [vcf_file] [outfile_prefix] [chromosome] [start] [end]
#     vcf_file: the file to convert to LDhat input format
#     outfile_prefix: prefix for all output files
#     chromosome: The name of the chromosome in the VCF file to process
#     start: start bp (inclusive) of region to make file for
#     end: end bp (inclusive) of region to make file for

from sys import argv
import gzip

infile_name = argv[1]
outfile_pref = argv[2]
sel_chr = argv[3]
reg_start = int(argv[4])
reg_end = int(argv[5])

reg_size = float(reg_end - reg_start)

infile = gzip.open(infile_name, "rt")
var_pos = []
n_snps = 0
for line in infile:
    # skip headers
    if line.startswith("##"):
        pass
    # process sample names
    if line.startswith("#CHROM"):
        header_fields = line.split()
        sample_start_ind = header_fields.index("FORMAT") + 1
        sample_end_ind = len(header_fields)
        sample_list = header_fields[sample_start_ind:sample_end_ind]
        gt_dict = {k: [] for k in sample_list} # initialize a dictionary of lists, sample names as keys
    # process the rest of the lines
    # for each line, loop through each sample
    
    if line.startswith(sel_chr):
        line_fields = line.split()
        # check that entry is within region specified
        if (int(line_fields[1]) >= reg_start and int(line_fields[1]) <= reg_end):
            pos_bp = int(line_fields[1]) # in bp
            rel_pos_bp = pos_bp - reg_start 
            var_pos.append(str(rel_pos_bp))
            n_snps += 1
            for i in range(sample_start_ind, sample_end_ind):
                # retrieve relevant fields from line for individual sample
                sample = sample_list[(i - sample_start_ind)]
                record = line.split()[i]
                record_gt = record.split(":")[0]
                # parse genotype into LDhat numeric code
                if record_gt == "0/0":
                    out_code = "0"
                elif record_gt == "1/1":
                    out_code = "1"
                elif record_gt in ["0/1", "1/0"]:
                    out_code = "2"
                else:
                    out_code = "?"
                # add to dict under correct sample
                gt_dict[sample].append(out_code)
infile.close()
# create output files
# .sites file
sites_outfile = open("{PREFIX}_{CHR}_{START}-{END}.sites".format(PREFIX = outfile_pref, CHR = sel_chr, START = reg_start, END = reg_end), "w")
sites_header = "{NSAMPLES}\t{NSNPS}\t{PHASED}\n".format(NSAMPLES = len(gt_dict.keys()), NSNPS = n_snps, PHASED = 2)
sites_outfile.write(sites_header)
for key in gt_dict.keys():
    gt_string = "".join(gt_dict[key])
    # insert linebreaks
    h = 0
    i = 80
    keep_going = True
    while keep_going:
        gt_string = gt_string[:i] + "\n" + gt_string[i:]
        h = i + 1
        i = h + 80
        keep_going = i < len(gt_string)
    sites_outfile.write(">{SAMPLE}\n".format(SAMPLE = key))
    sites_outfile.write(gt_string + "\n")
sites_outfile.close()
locs_outfile = open("{PREFIX}_{CHR}_{START}-{END}.locs".format(PREFIX = outfile_pref, CHR = sel_chr, START = reg_start, END = reg_end), "w")
# .locs file
locs_header = "{NSNPS}\t{CHR_SIZE}\t{MODEL}\n".format(NSNPS = str(n_snps), CHR_SIZE = reg_size, MODEL = "L")
locs_outfile.write(locs_header)
locs_outfile.write("\n".join(var_pos))
locs_outfile.close()
