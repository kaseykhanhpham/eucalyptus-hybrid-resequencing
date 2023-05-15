# script to convert genotype counts from VCFTools --counts tool to Ancestry_HMM input file format.
# assumes that the genetic map is already sorted by chromosome and position
# assumes that count files all have the same order and position of variants (i.e., they all came from the same VCF file)
# assumes that count files are consistent in which allele they order first and second
# should be agnostic to the number of references and admixed samples

# Usage: python make_ahmm_in.py [ref_filelist] [ref_min_cov] [admix_filelist] [genetic_dists] [out_name]
# ref_filelist: a file with the address of one reference file (the output of vcftools --counts) per line
# ref_min_cov: a file with the minimum number of alleles for each reference (corresponding with the order of ref_filelist) for a variant to be included
# admix_filelist: a file with the address of one admixed sample file (the output of vcftools --counts) per line
# genetic_dists: Table with the chromosome, position (bp), and genetic distance from the previous variant in that columnar order, sorted by chromosome and position. Can be generated from a genetic map and a counts file using accompanying script get_gen_dists.r

import sys

# Get CLI arguments
ref_filename = sys.argv[1]
cov_filename = sys.argv[2]
admix_filename = sys.argv[3]
dist_name = sys.argv[4]
out_name = sys.argv[5]

## Process reference count files ##
###################################

# open reference count files
ref_conns = []
ref_file = open(ref_filename, "r")
for line in ref_file:
    ref_conns.append(open(line.strip(), "r"))

# read min coverage (in individuals) for each ref
min_cov = []
cov_file = open(cov_filename, "r")
for line in cov_file:
    min_cov.append(int(line.strip()))
cov_file.close()

# skip first line (header) of each file
for conn in ref_conns:
    out_catcher = conn.readline()

allele1_ref_count = {} # dictionary -- keys are tuples (Chr, pos) and vals are list [ref1, ref2, ...]
allele2_ref_count = {} # dictionary -- keys are tuples (Chr, pos) and vals are list [ref1, ref2, ...]

for ref1_line in ref_conns[0]:
    skip = False
    allele1_temp = []
    allele2_temp = []

    # parse line from ref 1
    ref1_fields = ref1_line.strip().split()
    count_key = (ref1_fields[0], int(ref1_fields[1]))
    coverage = int(ref1_fields[3])
    allele1_temp.append(int(ref1_fields[4].split(":")[1]))
    allele2_temp.append(int(ref1_fields[5].split(":")[1]))

    # check for proper coverage
    if coverage < min_cov[0]:
        skip = True
    
    # now iterate through the rest of the reference files doing the same
    for i in range(1, len(ref_conns)):
        ref_line = ref_conns[i].readline()
        ref_fields = ref_line.strip().split()
        coverage = int(ref_fields[3])
        allele1_temp.append(int(ref_fields[4].split(":")[1]))
        allele2_temp.append(int(ref_fields[5].split(":")[1]))

        if coverage < min_cov[i]:
            skip = True
    
    # record reference counts for the variant if all coverages were high enough
    if not skip:
        allele1_ref_count[count_key] = allele1_temp
        allele2_ref_count[count_key] = allele2_temp

# close reference count files
for conn in ref_conns:
    conn.close()

## Process admixed samples ##
#############################

# Open admixed sample files
admix_conns = []
admix_file = open(admix_filename, "r")
for line in admix_file:
    admix_conns.append(open(line.strip(), "r"))

# skip header line of admixed count files
for conn in admix_conns:
    out_catcher = conn.readline()

allele1_adm_count = {} # dictionary -- keys are tuples (Chr, pos) and vals are list [ref1, ref2, ...]
allele2_adm_count = {} # dictionary -- keys are tuples (Chr, pos) and vals are list [ref1, ref2, ...]

for adm1_line in admix_conns[0]:
    allele1_temp = []
    allele2_temp = []

    # store information for first sample file
    adm1_fields = adm1_line.strip().split()
    count_key = (adm1_fields[0], int(adm1_fields[1]))
    allele1_temp.append(int(adm1_fields[4].split(":")[1]))
    allele2_temp.append(int(adm1_fields[5].split(":")[1]))

    # do the same for the rest of the admixed samples
    for i in range(1, len(admix_conns)):
        adm_line = admix_conns[i].readline()
        adm_fields = adm_line.strip().split()
        allele1_temp.append(int(adm_fields[4].split(":")[1]))
        allele2_temp.append(int(adm_fields[5].split(":")[1]))
    
    # record counts only if the variant passed reference coverage filters
    if count_key in allele1_ref_count.keys():
        allele1_adm_count[count_key] = allele1_temp
        allele2_adm_count[count_key] = allele2_temp

# Close admixed sample files
for conn in admix_conns:
    conn.close()

## Get genetic distances ##
###########################

dist_conn = open(dist_name, "r")

# Retrieve distances between remaining variants
gen_dists = []

dist_counter = 0
for line in dist_conn:
    chr = line.strip().split()[0]
    pos = int(line.strip().split()[1])
    dist = float(line.strip().split()[2])
    dist_counter = dist_counter + dist

    # check if variant was included in final set
    if (chr, pos) in allele1_ref_count.keys():
        # if so, record distance to final list and reset dist_counter for next variant
        gen_dists.append(dist_counter)
        dist_counter = 0

dist_conn.close()

## Construct final table ##
###########################

sorted_pos = sorted(allele1_ref_count.keys())
outfile = open(out_name, "w")

for i in range(0,len(sorted_pos)):
    out_line = sorted_pos[i][0] + "\t" + str(sorted_pos[i][1]) + "\t"
    # iterate over reference panels and add allele counts for each
    for j in range(0,len(ref_conns)):
        out_line = out_line + str(allele1_ref_count[sorted_pos[i]][j]) + "\t"
        out_line = out_line + str(allele2_ref_count[sorted_pos[i]][j]) + "\t"
    # add genetic distances between variants
    out_line = out_line + str(round(gen_dists[i], 3)) + "\t"
    # iterate over admixed samples and add allele counts for each
    for k in range(0, len(admix_conns)):
        out_line = out_line + str(allele1_adm_count[sorted_pos[i]][k]) + "\t"
        out_line = out_line + str(allele2_adm_count[sorted_pos[i]][k])
        # only add tab if not end of the line
        if k < len(admix_conns) - 1:
            out_line = out_line + "\t"
        else:
            out_line = out_line + "\n"
    # write to output file
    outfile.write(out_line)

outfile.close()
