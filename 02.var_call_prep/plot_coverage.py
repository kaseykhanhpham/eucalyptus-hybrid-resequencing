#!/usr/bin/env python

# python script to plot coverage of whole genome shotgun sequencing
# usage: python plot_coverage.py [infile_name] [sample_name] [average_outfile_name]

import sys
import numpy as np
from matplotlib import pyplot as plt

infile_name = sys.argv[1]
sample_name = sys.argv[2]
avg_outfile_name = sys.argv[3]

# initialize storage variables
cvg_sum = 0
num_rows = 0
sw_pos_sum = 0
sw_cvg_sum = 0
sw_pos = []
sw_cvg = []
infile = open(infile_name, "r")

# process input depth file
for line in infile:
    # parse line
    pos = int(line.split()[1])
    cvg = int(line.split()[2]) + int(line.split()[3])
    # store per-sample values
    cvg_sum = cvg_sum + cvg
    num_rows = num_rows + 1
    # store sliding window values
    sw_pos_sum = sw_pos_sum + pos
    sw_cvg_sum = sw_cvg_sum + cvg
    # every 100,000 rows, store average of values over that interval
    if num_rows % 100000 == 0:
        # calculate averages
        sw_pos_avg = sw_pos_sum / 100000
        sw_cvg_avg = sw_cvg_sum / 100000
        # store averages
        sw_pos.append(sw_pos_avg)
        sw_cvg.append(sw_cvg_avg)
        # reset sliding window variables
        sw_pos_sum = 0
        sw_cvg_sum = 0

infile.close()

# write averages to output averages file
avg_cvg = cvg_sum / num_rows
avg_outfile = open(avg_outfile_name, "a")
avg_outfile.write("{sample}\t{coverage}\n".format(sample = sample_name, coverage = avg_cvg))
avg_outfile.close()

# plot and export sliding window coverage
xvals = np.array(sw_pos)
yvals = np.array(list(map(np.log10, sw_cvg)))

plt.plot(xvals, yvals)
plt.savefig("{sample}_cover.png".format(sample = sample_name))
