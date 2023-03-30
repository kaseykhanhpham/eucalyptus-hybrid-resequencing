# Genome Scan Analysis

## Calculate nucleotide diversity statistics
First, I manually created lists of samples from the overall variant set that belonged to each population, under the names `Ecordata.txt`, `Eglobulus_MR.txt`, and `Eglobulus_ref.txt`. Each contains the RAPiD IDs under which each sample is identified in the VCF file, one per line, for each population. These will be used in all subsequent genome scan analyses to delimit populations.

Used [`vcftools`](https://vcftools.github.io) to calculate mean pairwise distance (pi) and Tajima's D in windows across the genome.

## Fst Calculations

**Calculate Fst in sliding windows along the genome:**

Calculated for Meehan Range _E. globulus_ and reference _E. globulus_ because I wanted to detect divergent regions in Meehan Range _E. globulus_ that may resemble _E. cordata_.

```bash
# Performed in UFRC queue system. See fst_windows.job for more details.
# Resources: 1.25 Mb, 1 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

vcftools --vcf "$INFILE" --out all_maf0.05 --weir-fst-pop "$POPLIST_DIR"/Eglobulus_MR.txt --weir-fst-pop "$POPLIST_DIR"/Eglobulus_ref.txt --fst-window-size 20000 --fst-window-step 2000
```

## Sliding Window Pi, Tajima's D, Dxy Calculations

Split multi-allelic SNPs into biallelic SNPs, remove problematic FORMAT fields which will cause errors:

```bash
# Performed in UFRC queue system. See format_egg.job for more details.
# Resources: 13 Mb, 2 min
module load bcftools/1.15

VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"
NAME="meehan_all_fil_maf0.00_snps"

bcftools annotate -x FORMAT/PL,FORMAT/AD,FORMAT/AO,FORMAT/QA -O v -o - "$VCFDIR"/"$NAME".vcf | bcftools norm -m -any -O v -o - > "$NAME"_biallelic.vcf
```

**Created list of chromosomes to iterate through**
Had to remove contigs which didn't have variants on them by hand before using this chromosome list file in the following script.

**Run python script to calculate Dxy in sliding windows across genome**
```bash
# Performed in UFRC queue system. See window_stats.job for more details.
# Resources: 210 Mb, 8 min
module load conda

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
NAME="meehan_all_fil_maf0.00_snps_biallelic"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/stats_windows.py "$WDIR"/"$NAME".vcf 100000 20000 "$WDIR"/pop_structure.json "$WDIR"/outgroup_structure.json "$WDIR"/chr_list.txt
```

## Patterson's D (and associated)

Use biallelic SNP set to calculate D statistics for individuals and in a sliding window across genomes, using _E. grandis_ as the outgroup.

```bash
# Performed in UFRC queue system. See dsuite.job for more details.
# Resources:
module load gcc/9.3.0
module load dsuite/0.4-r49

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"
VCFDIR="WDIR=/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
NAME="meehan_all_fil_maf0.00_snps_biallelic"

Dsuite Dinvestigate "$VCFDIR"/"$NAME".vcf "$WDIR"/SETS.txt "$WDIR"/test_trios.txt
```

## Analyze results
Get genome-wide distribution of FST and outlier regions (+2 sd).

```R
# Performed on local computer.
library(dplyr)

# Import files
wd <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/"
setwd(wd)

fst_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/windows_output/no_outgroup/all_maf0.00_no_outgr.windowed.weir.fst"
egglib_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/windows_output/no_outgroup/meehan_all_fil_maf0.00_snps_biallelic_formatted.vcf_scan_stats.tab"

fst_tab <- read.table(fst_filename, header = TRUE)
egglib_tab <- read.table(egglib_filename, header = TRUE, sep = "\t")

# Explore FST distribution
summary(fst_tab$MEAN_FST) # median = 0.005202, mean = 0.012169
sd(fst_tab$MEAN_FST) # sd = 0.02599
hist(fst_tab$MEAN_FST) # strong right skew
flagged_fst <- fst_tab$MEAN_FST > (mean(fst_tab$MEAN_FST) + 2*sd(fst_tab$MEAN_FST))
flagged_fst_windows <- as.data.frame(cbind("chr" = fst_tab[flagged_fst, "CHROM"], "start" = fst_tab[flagged_fst, "BIN_START"] - 1, "end" = fst_tab[flagged_fst, "BIN_END"])) # format to match egglib bin delimitation
flagged_fst_windows$start <- as.numeric(flagged_fst_windows$start)
flagged_fst_windows$end <- as.numeric(flagged_fst_windows$end)
# 1,978 windows of 54,541 were FST outliers

# Explore pi distribution
summary(egglib_tab$Pi) # median = 107.39, mean = 113.49
sd(egglib_tab$Pi, na.rm = TRUE) # sd = 99.64854
hist(egglib_tab$Pi) # has a bunch of small values which are probably pulling mean down a bit, but mean is still larger than median, so right skew is still the more dominant pattern..?
flagged_pi <- which(egglib_tab$Pi > (mean(egglib_tab$Pi, na.rm = TRUE) + 2*sd(egglib_tab$Pi, na.rm = TRUE)))
flagged_pi_windows <- egglib_tab[flagged_pi, c("chr", "start", "end")]
# 1,732 windows of 60,245 windows (54,703 w/o NAs) were pi outliers

# Explore dxy distribution
summary(egglib_tab$Dxy) # median = 0.220, mean = 0.226
sd(egglib_tab$Dxy, na.rm = TRUE) # sd = 0.05930694
hist(egglib_tab$Dxy) # more normal-looking than the others so far! Slight right skew.
flagged_dxy <- which(egglib_tab$Dxy > (mean(egglib_tab$Dxy, na.rm = TRUE) + 2*sd(egglib_tab$Dxy, na.rm = TRUE)))
flagged_dxy_windows <- egglib_tab[flagged_dxy, c("chr", "start", "end")]
# 1,908 windows of 60,245 windows (54,703 w/o NAs) were dxy outliers

# Explore Tajima's D distribution
summary(egglib_tab$Deta) # median = -0.003, mean = 0.029
sd(egglib_tab$Deta, na.rm = TRUE) # sd = 0.6164168
hist(egglib_tab$Deta) # approximately normal
# for directional selection
flagged_deta_neg <- which(egglib_tab$Deta < (mean(egglib_tab$Deta, na.rm = TRUE) - 2*sd(egglib_tab$Deta, na.rm = TRUE)))
flagged_deta_neg_windows <- egglib_tab[flagged_deta_neg, c("chr", "start", "end")]
# 1,058 windows of 60,245 windows (54,515 w/o NAs) were negative Tajima's D outliers
# for balancing selection
flagged_deta_pos <- which(egglib_tab$Deta > (mean(egglib_tab$Deta, na.rm = TRUE) + 2*sd(egglib_tab$Deta, na.rm = TRUE)))
flagged_deta_pos_windows <- egglib_tab[flagged_deta_pos, c("chr", "start", "end")]
# 1,946 windows of 60,245 windows (54,515 w/o NAs) were positive Tajima's D outliers

# Find overlaps in windows between outlier sets
fst_dxy_windows <- intersect(flagged_fst_windows, flagged_dxy_windows)
pi_dxy_windows <- intersect(flagged_pi_windows, flagged_dxy_windows)
dxy_deta_windows <- intersect(flagged_dxy_windows, flagged_deta_neg_windows)
```