# Python script to check for duplicate calls of the same position in a VCF file. Assumes VCF file is sorted by base pair position.
# Usage: python check_vcf_duplicates.py [infile] [out first] [out duplicates] [VCF without dupes]
#     infile: the VCF file to be checked
#     out first: file name for all of the first calls for a position
#     out duplicates: file name for all of the duplicate calls for a position
#     VCF without dupes: name for VCF output with duplicates removed

from sys import argv

infile_name = argv[1]
outfirst_name = argv[2]
outdup_name = argv[3]
vcf_out_name = argv[4]

infile = open(infile_name, "r")
outfirst = open(outfirst_name, "w")
outdup = open(outdup_name, "w")
vcf_out = open(vcf_out_name, "w")

# initiate storage variables for comparison
prev_line = ""
prev_chr = ""
prev_pos = 0

for line in infile:
    if line.strip().startswith("#"):
        # write header lines to output file automatically
        vcf_out.write(line)
    else:
        # splitting on whitespace, extract the first two columns
        var_chr = line.strip().split()[0]
        var_pos = line.strip().split()[1]

        # compare with previous line for duplicate position, indicating a duplicate call
        if var_chr == prev_chr and var_pos == prev_pos:
            # write to record output files
            outfirst.write(prev_line)
            outdup.write(line)
        else:
            # if not a duplicated line, write to output VCF
            vcf_out.write(line)

        # update storage variables for next iteration
        prev_line = line
        prev_chr = var_chr
        prev_pos = var_pos

infile.close()
outfirst.close()
outdup.close()
vcf_out.close()