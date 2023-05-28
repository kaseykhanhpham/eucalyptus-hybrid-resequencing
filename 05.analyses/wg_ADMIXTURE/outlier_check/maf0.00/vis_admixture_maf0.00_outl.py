###############
# Post-processing to visualize results of multiple ADMIXTURE runs and K values
# Run locally on personal computer because pong has interactive web visualization
# conda env: euc_hyb_reseq (python3.10, pip, scipy, numpy, plotnine, pong)

import pandas as pd

file_prefix = "all_fil_maf0.00_outl"

# Create pong filemap file with each value of K and each run
K_max = 6
n_runs = 10
qfile_dir = "C:\\Users\\Kasey\\OneDrive - University of Florida\\Grad School Documents\\Projects\\eucalyptus-hybrid-resequencing\\05.analyses\\wg_admixture\\outlier_check\\maf0.00\\admixture_output"

outfile = open("{PREFIX}_filemap.txt".format(PREFIX=file_prefix), "w")

for k in range(2,(K_max + 1)):
    for r in range(1,(n_runs + 1)):
        outfile.write("k{K}r{R}\t{K}\t{INDIR}\\{PREFIX}.K{K}.r{R}.Q\n".format(K=k, R=r, INDIR=qfile_dir, PREFIX=file_prefix))

outfile.close()

# Create pong ind2pop file for samples
spp_tab_loc = "C:\\Users\\Kasey\\OneDrive - University of Florida\\Grad School Documents\\Projects\\eucalyptus-hybrid-resequencing\\00.metadata\\03.seq_analysis\\sample_spp_table.csv"
fam_file_loc = "C:\\Users\\Kasey\\OneDrive - University of Florida\\Grad School Documents\\Projects\\eucalyptus-hybrid-resequencing\\05.analyses\\wg_admixture\\outlier_check\\maf0.00\\" + file_prefix + ".fam"

spp_table = pd.read_csv(spp_tab_loc, header = 0, index_col="RAPiD_ID")
spp_order = list(pd.read_table(fam_file_loc, sep = " +", header = None, engine = "python").iloc[:,0])
spp_table.loc[spp_order].to_csv("{PREFIX}_ind2pop.txt".format(PREFIX=file_prefix), sep = "\t", index = False, header = False, columns=["Taxon"])

# Create pong poporder file
outfile = open("{PREFIX}_poporder.txt".format(PREFIX=file_prefix), "w")
outfile.write("cord_MR\nglob_MR\nglob_pure\n")
outfile.close()

# Create pong color file
outfile = open("{PREFIX}_colors.txt".format(PREFIX=file_prefix), "w")
outfile.write("#FFA716\n#1A688A\n#E6E2DF\n#AF1E36\n#68C7D8\n#272623\n")
outfile.close()
