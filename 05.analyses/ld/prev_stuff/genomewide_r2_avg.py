# python script to average across all chromosomes worth of average r2 values per distance between base pairs
# only meant for this one instance of use in ld calculations
# usage: python genomewide_r2_avg.py [FILENAME_PATTERN]
#        where FILENAME_PATTERN is formatted so "{CHR}" goes where the chromosome name is supposed to

import pandas as pd
from sys import argv

pattern = argv[1]
chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11"]

conn_list = []

# open a file connection for each chromosome (corresponds in position to chr_list)
for chr in chr_list:
    conn_list.append(open(pattern.format(CHR = chr), "r"))

# skip header (first line in each file)
for conn in conn_list:
    conn.readline()

genomewide_r2 = pd.DataFrame(data = None, columns = ["dist", "r2"])

# read r2 for each distance for each chromosome simultaneously and average
for chr01_line in conn_list[0]:
    r2_list = []

    # process opened chr01 line first
    chr01_r2 = float(chr01_line.strip().split(",")[1])
    # check for stand-in NA value
    if chr01_r2 != -1:
        r2_list.append(chr01_r2)
    # process the rest of the chromosomes in the connection list
    for conn in conn_list[1:]:
        parsed_r2 = float(conn.readline().strip().split(",")[1])
        # check for stand-in NA value
        if parsed_r2 != -1:
            r2_list.append(parsed_r2)

    dist = chr01_line.strip().split(",")[0]
    avg_r2 = sum(r2_list)/len(r2_list) # take average for however many chr r2 values could be retrieved

    genomewide_r2.loc[len(genomewide_r2.index)] = [dist, avg_r2]

# close connections
for conn in conn_list:
    conn.close()

genomewide_r2.to_csv("genomewide_r2.csv", index = False)