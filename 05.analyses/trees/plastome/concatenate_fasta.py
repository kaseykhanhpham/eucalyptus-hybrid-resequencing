# python script to concatenate multiple FASTA files by sequence name
# usage: python concatenate_fasta.py [comma-delimited list of input FASTA files] [output name]
# WARNING: cannot handle FASTA files with comments

import sys

infile_list_str = sys.argv[1]
outfile_name = sys.argv[2]

# PARSE AND READ INPUT FILES
infile_name_list = infile_list_str.split(",")
infile_list = []
infile_len = []

for name in infile_name_list:
    infile_dict = {} # save sequences to a dictionary where keys are seq names
    infile_conn = open(name, "r")
    # read first header and initialize storage variables
    seqname = infile_conn.readline().strip(">").strip()
    seq = ""
    # loop through rest of file to read in seqs
    for line in infile_conn:
        if line.startswith(">"):
            # save previous sequence to dictionary
            infile_dict[seqname] = seq
            # reset storage variables
            seqname = line.strip(">").strip()
            seq = ""
        else:
            # append sequence line to storage variable
            seq = seq + line.strip()
    # write final sequence to dictionary
    infile_dict[seqname] = seq
    infile_conn.close()
    # save length of final sequence
    infile_len.append(len(seq))
    # append dictionary of seqs to storage list
    infile_list.append(infile_dict)

# APPEND SEQS BY NAME
# get unique list of sequences to append
seqname_list = []
for infile in infile_list:
    seqname_list = seqname_list + list(infile.keys())
seqname_set = set(seqname_list)

# iterate through unique sequences and concatenate sequences
concatenated_seqs_dict = {}
for seqname in seqname_set:
    seq_list = []
    # read in sequences from each input file
    for i in range(len(infile_list)):
        try:
            # extract sequence from the target file
            seq_list.append(infile_list[i][seqname])
        except KeyError:
            # add null dashes if sequence name doesn't exist for the target file
            insert_seq = "-" * infile_len[i]
            seq_list.append(insert_seq)
    # concatenate all extracted sequences for the sequence name
    concat_seq = ""
    for seq in seq_list:
        concat_seq = concat_seq + seq
    # save the concatenated sequence to dictionary
    concatenated_seqs_dict[seqname] = concat_seq

# WRITE FASTA FILE
outfile_conn = open(outfile_name, "w")
for seqname in concatenated_seqs_dict:
    outfile_conn.write(">{SEQNAME}\n".format(SEQNAME = seqname))
    outfile_conn.write("{SEQ}\n".format(SEQ = concatenated_seqs_dict[seqname]))
outfile_conn.close()