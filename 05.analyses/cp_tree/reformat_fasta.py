# -*- coding: utf-8 -*-
"""
Function: Split FASTA sequences on one line to multiple lines
Author: Kasey Pham
Date edited: 05/31/2020
Usage: python reformat_fasta.py [num characters per line] [input name] [output name]
"""

import sys

numchar = int(sys.argv[1])
infile_name = sys.argv[2]
outfile_name = sys.argv[3]

infile = open(infile_name, "r")
outfile = open(outfile_name, "w")
counter = 0
while True:
    char_stream = infile.read(1)
    if not char_stream:
        print("End of file reached.")
        break
    if char_stream == ">":
        outfile.write(">")
        outfile.write(infile.readline())
        counter = 0
    else:
        if counter == numchar:
            outfile.write("\n")
            counter = 0
        outfile.write(char_stream)
        counter = counter + 1
infile.close()
outfile.close()