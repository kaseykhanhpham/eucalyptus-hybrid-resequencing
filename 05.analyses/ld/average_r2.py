# python script to average r2 for each pairwise distance between SNPs in an emeraLD r2 output file
# most appropriate to run via job scheduler since it will require about as much RAM as the size of the emeraLD output file.

from statistics import mean
import sys
import pandas as pd
import gzip
import argparse

# parse user-passed arguments
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", help = "LD calculation file from emeraLD", type = str)
parser.add_argument("-o", "--outfile", help = "Name of output average table", type = str)

args = parser.parse_args()

infile_name_stripped = args.infile.split("/")[(len(args.infile.split("/")) - 1)]
outfile_name_stripped = args.outfile.split("/")[(len(args.outfile.split("/")) - 1)]

# initialize storage dictionary
all_r2_dict = {}

# open input file
infile = gzip.open(args.infile, "rb")
# read header
infile.readline()

# loop through each line of the input file
for line in infile:
    line_cols = line.split()
    # get distance between SNPs (BP_B - BP_A)
    bp_a = int(line_cols[1])
    bp_b = int(line_cols[2])
    snp_dist = bp_b - bp_a

    # check if key exists for distance already in dictionary
    if not (snp_dist in all_r2_dict.keys()):
        # if it doesn't, create key with empty list as value
        all_r2_dict[snp_dist] = []
    # append R2 to list
    all_r2_dict[snp_dist].append(float(line_cols[4]))

# close infile
infile.close()

if len(all_r2_dict.keys()) == 0:
    print("{INFILE} is empty".format(INFILE = infile_name_stripped))
else:
    # initialize dictionary for storing mean r2s
    avg_r2_dict = {}
    # loop through dictionary of lists
    for snp_dist in range(1, (max(all_r2_dict.keys()) + 1)):
        try:
            # average r2 values and save to new dict
            avg_r2_dict[snp_dist] = mean(all_r2_dict[snp_dist])
        except KeyError:
            avg_r2_dict[snp_dist] = -1 # represent NAs as an impossible value for processing later

    # convert averaged dictionary to pandas dataframe
    avg_r2_df = pd.DataFrame(list(avg_r2_dict.items()), columns = ["dist", "r2"])

    # sort dataframe by SNP distance
    avg_r2_df.sort_values(by=["dist"], inplace = True)

    # export pandas dataframe
    avg_r2_df.to_csv(args.outfile, index = False)
    print("{INFILE} averages written to {OUTFILE}.".format(INFILE = infile_name_stripped, OUTFILE = outfile_name_stripped))