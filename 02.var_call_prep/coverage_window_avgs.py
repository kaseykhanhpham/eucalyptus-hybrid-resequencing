#!/usr/bin/env python

# python script to calculate average coverage of whole genome shotgun sequencing in sliding windows
# using the output of samtools depth
# usage: python coverage_window_avgs.py -i [infile] -o [outfile] -a [summaryfile] -s [windowsize] -n [sampleid]
#        infile: name of depth file from samtools to average
#        outfile: name of output file to write with sliding window averages for this sample
#        summaryfile: name of file in which to append overall average coverage for this specific sample
#                     (under assumption this script will be run multiple times). Will create if it doesn't
#                     already exist.
#        windowsize: sliding window size across which to average coverage
#        sampleid: ID of sample coverage is calculated for

import argparse
import subprocess

# parse user-passed arguments
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", help = "Coverage input file name", type = str)
parser.add_argument("-o", "--outfile", help = "Sliding window coverage output file name", type = str)
parser.add_argument("-a", "--summaryfile", help = "Multi-sample summary output file name", type = str, default = "average_coverage.txt")
parser.add_argument("-s", "--windowsize", help = "Averaging sliding window size", type = int, default = 100000)
parser.add_argument("-n", "--sampleid", help = "Name or ID of sample", type = str, default = "Sample")

args = parser.parse_args()

print("Calculating sliding window averages for sample {SAMPLE} from {INFILE} with window size {WINSIZE}.".format(SAMPLE = args.sampleid, INFILE = args.infile, WINSIZE = args.windowsize))

# initialize storage variables
cvg_sum = 0
num_rows = 0
sw_pos_sum = 0
sw_cvg_sum = 0
sw_pos = []
sw_cvg = []
infile = open(args.infile, "r")

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
    if num_rows % args.windowsize == 0:
        # calculate averages
        sw_pos_avg = sw_pos_sum / args.windowsize
        sw_cvg_avg = sw_cvg_sum / args.windowsize
        # store averages
        sw_pos.append(sw_pos_avg)
        sw_cvg.append(sw_cvg_avg)
        # reset sliding window variables
        sw_pos_sum = 0
        sw_cvg_sum = 0

infile.close()

# MULTI-SAMPLE SUMMARY FILE
avg_cvg = cvg_sum / num_rows
# Create summary file if it does not already exist at the specified location
subprocess.run(["touch", args.summaryfile])
# Write genome-wide average to summary file
print("Writing genome-wide average coverage to {SUMFILE}.".format(SUMFILE = args.summaryfile))
summary_outfile = open(args.summaryfile, "a")
summary_outfile.write("{SAMPLE}\t{COV:2f}\n".format(SAMPLE = args.sampleid, COV = avg_cvg))
summary_outfile.close()

# SAMPLE-SPECIFIC SLIDING WINDOW FILE
# Write sliding window averages to output file
print("Writing sliding window averages to {OUTFILE}.".format(OUTFILE = args.outfile))
outfile = open(args.outfile, "w")
for i in range(len(sw_pos)):
    outfile.write("{POS:0f}\t{COV:2f}\n".format(POS = sw_pos[i], COV = sw_cvg[i]))
outfile.close()
# 
# plot and export sliding window coverage
# xvals = np.array(sw_pos)
# yvals = np.array(list(map(np.log10, sw_cvg)))

# plt.plot(xvals, yvals)
# plt.savefig("{sample}_cover.png".format(sample = sample_name))
