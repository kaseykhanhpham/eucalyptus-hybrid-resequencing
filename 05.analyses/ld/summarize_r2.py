# python script to summarize linkage disequilibrium r2 tables after running average_r2.py
# Usage: python summarize_r2.py [infile name] [threshhold of r2]

import sys

infile_name = sys.argv[1]
threshhold = float(sys.argv[2])
infile = open(infile_name, "r")

# skip header
infile.readline()

# initialize counter variables
first = 0
last = 0
first_cross_status = False
curr_cross_status = False

# iterate through CSV file
for line in infile:
    dist = int(line.strip().split(",")[0])
    r2 = float(line.strip().split(",")[1])
    # record first time r2 crosses threshhold
    if not first_cross_status:
        if r2 <= threshhold:
            first_cross_status = True
            curr_cross_status = True
            first = dist
    else: 
        if r2 <= threshhold:
            # if threshhold is crossed + to - again after first time, update counters
            if not curr_cross_status:
                curr_cross_status = True
                last = dist
            else:
                # stop loop if r2 has been above threshhold for 500bp
                if last - dist > 499:
                    break
        # if threshhold is crossed + to - update counter
        else:
            curr_cross_status = False

# print summary
print("min LD at {THRESHHOLD}: {FIRST}".format(THRESHHOLD = threshhold, FIRST = first))
print("max LD at {THRESHHOLD}: {LAST}".format(THRESHHOLD = threshhold, LAST = last))
print("midpoint LD at {THRESHHOLD}: {MID}".format(THRESHHOLD = threshhold, MID = first + (last - first)/2))
