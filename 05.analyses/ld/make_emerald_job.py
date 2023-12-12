# python script to make emeraLD job which calculates LD for non-overlapping 100kb sliding windows across a chromosome
# usage: python make_emerald_job.py [vcf_loc] [chr_name] [chr_size] [window_size] [mac] [ind_file] [outfile_name]
#     vcf_loc: the address of the VCF to use in emeraLD
#     chr_name: ID of chromosome to generate commands for
#     chr_size: size of chromosome in base pairs to generate commands for
#     window_size: size of sliding window (non-overlapping) across the chromosome in base pairs
#     mac: minimum minor allele count
#     ind_file: file of individuals to include
#     outfile_name: name of job file to write

# Note: You will still need to edit the heading afterwards to be compatible with queue/module system of choice
# emeraLD has several dependencies, including python, tabix, and htslib

from sys import argv
from math import floor

# import arguments
vcf_loc = argv[1]
chr_name = argv[2]
chr_size = int(argv[3])
window_size = int(argv[4])
mac = argv[5]
ind_file = argv[6]
out_name = argv[7]

# Populate starts and ends
window_starts = []
window_ends = []

for i in range(floor(chr_size / window_size)):
    window_starts.append(i * window_size)

for j in window_starts[1:]:
    window_ends.append(j - 1)
window_ends.append(chr_size)

# write to job
header_template = "#!/bin/bash\n#SBATCH --account=soltis\n#SBATCH --qos=soltis-b\n#SBATCH --job-name=glob{CHR}m{MAC}_ld\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=kasey.pham@ufl.edu\n#SBATCH --mem=100mb\n#SBATCH --time=1:00:00\n#SBATCH --cpus-per-task=1\n#SBATCH --nodes=1\n#SBATCH --output=glob_windows_{CHR}_mac{MAC}_ld_%j.out\n#SBATCH --error=glob_windows_{CHR}_mac{MAC}_ld_%j.err\nmodule purge\nmodule load htslib/1.15\nmodule load emerald/0.1\n\n"
command_template = "emeraLD --in {VCF} --phase --region {CHR}:{START}-{END} --mac {MAC} --include {INDS} --stdout | bgzip -c > {CHR}_{START}-{END}_ld.txt.gz\n"

out_conn = open(out_name, "w")
out_conn.write(header_template.format(CHR = chr_name, MAC = mac))
for k in range(len(window_starts)):
    out_conn.write(command_template.format(VCF = vcf_loc, CHR = chr_name, START = window_starts[k], END = window_ends[k], MAC = mac, INDS = ind_file))

out_conn.close()
