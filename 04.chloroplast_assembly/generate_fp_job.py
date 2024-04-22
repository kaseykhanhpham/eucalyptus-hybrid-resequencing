# python script to generate HiperGator job for FastPlast
# everything is hardcoded in since it's specific to my local specs

import pandas as pd

metadata_loc = "/blue/soltis/kasey.pham/euc_hyb_reseq/sample_sequencing_metadata.csv"
output_file = "/blue/soltis/kasey.pham/euc_hyb_reseq/cp_assembly/fastplast.job"
reads_dir = "/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
ref_loc = "/blue/soltis/kasey.pham/euc_hyb_reseq/refs/organelle/AY780259.1"

meta_tab = pd.read_csv(metadata_loc)

# Create dictionary where key is Sample ID and values are lists of RunSample IDs
# Mask to just R1
dir_mask = meta_tab["Direction"] == "R1"
# Get unique list of Sample IDs
sample_ids = meta_tab["Sample"].unique()
# Loop through Sample IDs and mask to get corresponding RunSample IDs
sample_dic = {}
for sample in sample_ids:
    sample_mask = meta_tab[dir_mask]["Sample"] == sample
    sample_dic[sample] = list(meta_tab[dir_mask][sample_mask]["RunSample"])

# Write job submission header
job_text = '''#!/bin/bash
#SBATCH --account=soltis
#SBATCH --qos=soltis-b
#SBATCH --job-name=fastplast
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kasey.pham@ufl.edu
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1

pwd; hostname; date
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

module purge
module load fastplast/1.2.8

'''

# Iterate through dictionary and write values to job body
for sample in sample_dic:
    fp_line = "perl $HPC_FASTPLAST_DIR/Fast-Plast/fast-plast.pl -1 {READS_DIR}/{RUNSAMPLE1}_R1_paired.fq,{READS_DIR}/{RUNSAMPLE2}_R1_paired.fq -2 {READS_DIR}/{RUNSAMPLE1}_R2_paired.fq,{READS_DIR}/{RUNSAMPLE2}_R2_paired.fq -n {SAMPLE} --subsample 30000000 --threads 7 --user_bowtie {REF_LOC} --clean deep --skip trim --coverage_analysis --min_coverage 10 \n\n".format(READS_DIR = reads_dir, RUNSAMPLE1 = sample_dic[sample][0], RUNSAMPLE2 = sample_dic[sample][1], SAMPLE = sample, REF_LOC = ref_loc)
    job_text = job_text + fp_line

# Write job
outp_conn = open(output_file, "w")
outp_conn.write(job_text)
outp_conn.close()
