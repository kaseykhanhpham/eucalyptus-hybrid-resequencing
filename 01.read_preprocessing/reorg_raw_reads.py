#!/usr/bin/env python

# python script to create symlinks of all raw reads in new directory
# Kasey Pham, Jan 11 2022, Eucalyptus globulus x cordata resequencing project

import os

indir = "/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/RAPiD_raw_reads"
outdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/reads/raw_reads"
meta_tbl_name = "/blue/soltis/kasey.pham/euc_hyb_reseq/sample_sequencing_metadata.csv"
inname_col = 0
outname_col = 8
dir_col = 10

meta_tbl = open(meta_tbl_name, "r")

# skip header line
meta_tbl.readline()

# iterate through rest of metadata table
for line in meta_tbl:
    # parse fields in line
    line_list = line.split(",")
    infile = line_list[inname_col]
    seqid = line_list[outname_col]
    direction = line_list[dir_col]
    # create symlink from parsed fields
    os.system("ln -s {INDIR}/{INFILE} {OUTDIR}/{SEQID}_{DIRECTION}_raw.fastq.gz".format(INDIR = indir, INFILE = infile, OUTDIR = outdir, SEQID = seqid, DIRECTION = direction))

meta_tbl.close()
