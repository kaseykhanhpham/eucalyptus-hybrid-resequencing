# -*- coding: utf-8 -*-
"""
Description: Python code to generate Picard commands to add read groups to SAM files
             given a table of sample names and corresponding read group from stdin.
             Outputs a job file of all Picard commands to add read group to each
             sample in the table.
Author: Kasey
Last Updated: 01/20/2022

Usage: python make_rg_job.py <readgroup table> <input directory> <input extension> <output directory> <output extension> <account_name> <email>
Arguments:
    job file name: the name of the job file to be written. Overrides any file at
                   that address.
"""

import sys

readgroup_tb = sys.argv[1]
input_dir = sys.argv[2]
input_ext = sys.argv[3]
output_dir = sys.argv[4]
output_ext = sys.argv[5]
account_name = sys.argv[6]
email = sys.argv[7]

sample_table = open(readgroup_tb, mode = "r")
job_file = open("add_readgroups.job", mode = "w")
header_text = '''#!/bin/bash
#SBATCH --account={ACCOUNT_NAME}
#SBATCH --qos={ACCOUNT_NAME}-b
#SBATCH --job-name=add_rg
#SBATCH --mail-type=ALL
#SBATCH --mail-user={EMAIL}
#SBATCH --mem=15gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=add_rg_%j.out
#SBATCH --error=add_rg_%j.err

module load picard/2.25.5

'''.format(ACCOUNT_NAME = account_name, EMAIL = email)

job_file.write(header_text)

# skip header
sample_table.readline()

for row in sample_table:
    # parse table row
    sample = row.split()[0]
    barcode = row.split()[1]
    run_sample = row.split()[2]
    flowcell = row.split()[5]
    lane = row.split()[6]
    
    command = "picard AddOrReplaceReadGroups -I {INPUT_DIR}/{RUN_SAMPLE}{INPUT_EXT} -O {OUTPUT_DIR}/{RUN_SAMPLE}{OUTPUT_EXT} -LB {SAMPLE}.leaf -PL ILLUMINA -PU {FLOWCELL}.{LANE}.{BARCODE} -SM {SAMPLE} -ID {FLOWCELL}.{LANE}.{BARCODE}\n".format(INPUT_DIR = input_dir, RUN_SAMPLE = run_sample, INPUT_EXT = input_ext, OUTPUT_DIR = output_dir, OUTPUT_EXT = output_ext, SAMPLE = sample, FLOWCELL = flowcell, LANE = lane, BARCODE = barcode)
    job_file.write(command)
sample_table.close()
job_file.close()
print("Done writing to add_readgroups.job!")