# python script to make emeraLD job which calculates LD for non-overlapping 100kb sliding windows across a chromosome
# usage: python make_emerald_job.py -v [vcf_loc] -c [chr_name] -s [chr_size] -i [ind_file] -o [outfile_name] -w [window_size] -m [mac]
#     vcf_loc: the address of the VCF to use in emeraLD
#     chr_name: ID of chromosome to generate commands for
#     chr_size: size of chromosome in base pairs to generate commands for
#     window_size: size of sliding window (non-overlapping) across the chromosome in base pairs
#     mac: minimum minor allele count
#     ind_file: file of individuals to include
#     outfile_name: name of job file to write

# Note: You will still need to edit the heading afterwards to be compatible with queue/module system of choice
# emeraLD has several dependencies, including python, tabix, and htslib

from math import floor
import argparse

# parse user-passed arguments
parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", help = "VCF input file", type = str)
parser.add_argument("-c", "--chr", help = "Name of focal chromosome", type = str)
parser.add_argument("-s", "--size", help = "Chromosome size", type = int)
parser.add_argument("-i", "--indfile", help = "File of individuals to include", type = str)
parser.add_argument("-o", "--outfile", help = "Name of output file", type = str)
parser.add_argument("-t", "--tax", help = "Name of the taxon", type = str, default = "tax")
parser.add_argument("-w", "--winsize", help = "Size of sliding windows", type = int, default = 100000)
parser.add_argument("-m", "--mac", help = "Minor allele frequency", type = int, default = 1)


args = parser.parse_args()

# Populate starts and ends
window_starts = []
window_ends = []

for i in range(floor(args.size / args.winsize)):
    window_starts.append(i * args.winsize)

for j in window_starts[1:]:
    window_ends.append(j - 1)
window_ends.append(args.size)

# write to job
header_template = "#!/bin/bash\n#SBATCH --account=soltis\n#SBATCH --qos=soltis-b\n#SBATCH --job-name={TAX}{CHR}_ld\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=kasey.pham@ufl.edu\n#SBATCH --mem=100mb\n#SBATCH --time=1:00:00\n#SBATCH --cpus-per-task=1\n#SBATCH --nodes=1\n#SBATCH --output={TAX}_windows_{CHR}_mac{MAC}_ld_%j.out\n#SBATCH --error={TAX}_windows_{CHR}_mac{MAC}_ld_%j.err\nmodule purge\nmodule load htslib/1.15\nmodule load emerald/0.1\n\n"
command_template = "emeraLD --in {VCF} --phase --region {CHR}:{START}-{END} --mac {MAC} --include {INDS} --stdout | bgzip -c > {CHR}_{START}-{END}_ld.txt.gz\n"

out_conn = open(args.outfile, "w")
out_conn.write(header_template.format(CHR = args.chr, MAC = args.mac, TAX = args.tax))
for k in range(len(window_starts)):
    out_conn.write(command_template.format(VCF = args.vcf, CHR = args.chr, START = window_starts[k], END = window_ends[k], MAC = args.mac, INDS = args.indfile))

out_conn.close()
