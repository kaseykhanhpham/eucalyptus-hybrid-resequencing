# dadi demography model fitting
# step: site frequency spectrum construction and bootstrapping
# relies on code from dadi manual examples: https://dadi.readthedocs.io/en/latest/examples/fs_from_data/fs_from_data/

# Pickle is used to save variables as files for future use
import pickle
# MatPlotLib is a libary dadi uses for plotting frequency spectrum
import matplotlib.pyplot as plt
import dadi
# import os

# names for output files
wdir = "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/dadi"
name_stem = "globMR_cordMR_ns32-16"

# Parse the VCF file to generate a data dictionary
datafile = "/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil_biallelic.prune.vcf"
dd = dadi.Misc.make_data_dict_vcf(datafile, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/demography/stairway_plot/esfs_poplist.txt")

# Extract the spectrum for MR E. globulus and MR E. cordata from that dictionary, with 
# E. globulus projected down to 32 of 40 and E. cordata projected down to 16 of 20.
# haplotype counts at which the number of segregating sites was maximized
pop_ids, ns = ["glob_MR", "cord_MR"], [32, 16]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized = False) # folded since ancestral alleles not known

# Make directory for saving data dictionaries
# if not os.path.exists("{WDIR}/data/data_dicts".format(WDIR = wdir)):
#  os.makedirs("{WDIR}/data/data_dicts".format(WDIR = wdir))

# Saved data dictionary to disk
pick = open("{WDIR}/data/data_dicts/{NAME}.bpkl".format(WDIR = wdir, NAME = name_stem), "wb")
pickle.dump(dd, pick, 2)
# Saved extracted spectrum
fs.to_file("{WDIR}/data/fs/{NAME}.fs".format(WDIR = wdir, NAME = name_stem))

# Generate 100 bootstrap datasets, by dividing the genome into 500kb chunks and
# resampling from those chunks.
Nboot, chunk_size = 100, 5e5
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, ns, polarized = False)
# Saved bootstraps
for i in range(len(boots)):
    boots[i].to_file("{WDIR}/data/fs/bootstraps/{NAME}.boot_{NUM}.fs".format(WDIR = wdir, NAME = name_stem, NUM = str(i)))

# Plotted site frequency spectrum
fig = plt.figure(982342)
fig.clear()
dadi.Plotting.plot_single_2d_sfs(fs)
fig.savefig("{WDIR}/data/fs/{NAME}.png".format(WDIR = wdir, NAME = name_stem))
plt.close(fig)