# python script to assemble FASTA files of each chloroplast region from Fastplast assemblies

import os
sample_list_file = "/blue/soltis/kasey.pham/euc_hyb_reseq/sample_ids.txt"
cp_dir = "/blue/soltis/kasey.pham/euc_hyb_reseq/cp_assembly"
outgroups = ["HM347959.1", "KC180790.1"]

# Read in sample list
sample_list_conn = open(sample_list_file, "r")
sample_list = []
for line in sample_list_conn:
    sample_list.append(line.strip())
sample_list_conn.close()

# Split each cp region into its own file and rename FASTA header to sample ID
for sample in sample_list:
    in_file = "{CP_DIR}/{SAMPLE}/Final_Assembly/{SAMPLE}_CP_pieces.fsa".format(CP_DIR = cp_dir, SAMPLE = sample)
    in_conn = open(in_file, "r")
    out_conn = open("dummy_file", "w")
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

# Add outgroups
for outgroup in outgroups:
    in_file = "{CP_DIR}/{OUTGROUP}/{OUTGROUP}_CP_pieces.fsa".format(CP_DIR = cp_dir, OUTGROUP = outgroup)
    in_conn = open(in_file, "r")
    out_conn = open("dummy_file", "w")
    for line in in_conn:
        if line.startswith(">"):
            # close previously open file
            out_conn.close()
            # read current chloroplast region
            region = line.strip(">").strip()
            # open relevant region file and write outgroup ID as header
            out_conn = open("{REGION}_temp.fas".format(REGION = region), "a")
            out_conn.write(">{OUTGROUP}\n".format(OUTGROUP = outgroup))
        else:
            # just write line to currently open file if part of sequence body
            out_conn.write(line)
    in_conn.close()
out_conn.close()


os.system("rm dummy_file")