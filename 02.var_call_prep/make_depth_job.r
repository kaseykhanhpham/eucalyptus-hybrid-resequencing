# -*- coding: utf-8 -*-
# R script to generate job to calculate depth over multiple sequencing runs for samples
# given an input table with sample information and location of files
# Currently hardcoded with specific files, not generalizeable!
# Usage: Rscript make_depth_job.r

sample_table_name <- "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/05.depth/runid_to_sampleid.csv"
input_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
infile_suffix <- "_marked.bam"
sample_col <- 3
run_id_col <- 2

account_name <- "soltis"
email <- "kasey.pham@ufl.edu"

sample_table <- read.csv(sample_table_name, header = TRUE, as.is = TRUE)
out_sample_list <- unique(sample_table[,sample_col])

job_text <- "#!/bin/bash
#SBATCH --account=ACCOUNT_NAME
#SBATCH --qos=ACCOUNT_NAME-b
#SBATCH --job-name=depth
#SBATCH --mail-type=ALL
#SBATCH --mail-user=EMAIL$
#SBATCH --mem=10gb
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=depth_%j.out
#SBATCH --error=depth_%j.err

module load samtools/1.12

"

job_text <- gsub("ACCOUNT_NAME", account_name, job_text)
job_text <- gsub("EMAIL", email, job_text)

for(out_sample in out_sample_list){
    sample_ids <- sample_table[which(sample_table[,sample_col] == out_sample), run_id_col]
    command_text <- "samtools depth -o /blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/05.depth/NAME_cover.txt FILES\n"
    command_text <- gsub("NAME", out_sample, command_text)
    command_text <- gsub("FILES", paste(paste(input_dir, "/", sample_ids, infile_suffix, sep = ""), collapse = " "), command_text)
    job_text <- paste(job_text, command_text, "\n", sep = "")
}

write.table(job_text, "depth.job", quote = FALSE)