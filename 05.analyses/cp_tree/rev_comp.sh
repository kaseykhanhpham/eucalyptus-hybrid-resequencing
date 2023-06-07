#!/bin/bash

# sysinfo_page - a script to automate reverse-complementing entries in a FASTA file using `biopieces`
# usage: rev_comp.sh [name root] [entries to reverse-complement] [entries to ignore]
# name root: the name of the FASTA file to reverse-complement minus the ".fsa" file extension (currently only works with that extension)
# entries to reverse-complement: comma-separated
# entries to ignore: comma-separated
# note: assumes the `biopieces` program is installed and the module is loaded before running the script
# e.g., module load biopieces/2.0

#### Command line args
echo doing "$0".fsa

mkdir $BP_DATA $BP_TMP $BP_LOG

mv "$0".fsa "$NAME"_raw.fsa
# reverse complement just the specified entries
read_fasta -i "$0"_raw.fsa | grab -p "$1" | reverse_seq | complement_seq | write_fasta -x -o "$0".fsa
# extract the specified entries as-is
read_fasta -i "$0"_raw.fsa | grab -p "$2" | write_fasta -x -o temp_warningthisoverwritesthisfilename.fsa
# merge reverse-complemented sequences with the untouched sequences
cat temp_warningthisoverwritesthisfilename.fsa >> "$0".fsa
rm temp_warningthisoverwritesthisfilename.fsa