# R script to read in Liu et al.'s code for converting VCF to SFS and run on all three Eucalyptus populations.

# Read in code
script_dir <- "/blue/soltis/kasey.pham/bin"
w_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/stairway_plot"
source(paste(script_dir, "vcf2sfs.r", sep = "/"))

setwd(w_dir)

# Run for each population
# E. globulus introgressants
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "glob_MR.fs", 1)

# E. globulus reference
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "glob_pure.fs", 2)

# E. cordata reference
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "cord.fs", 3)
