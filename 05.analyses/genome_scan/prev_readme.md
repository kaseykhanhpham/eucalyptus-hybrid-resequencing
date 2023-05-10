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

vcftools --vcf "$INFILE" --out all_maf0.00 --weir-fst-pop "$POPLIST_DIR"/Eglobulus_MR.txt --weir-fst-pop "$POPLIST_DIR"/Eglobulus_ref.txt --fst-window-size 20000 --fst-window-step 2000
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

**Run python script to calculate stats in sliding windows across genome**
```bash
# Performed in UFRC queue system. See window_stats.job for more details.
# Resources: 210 Mb, 8 min
module load conda

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
NAME="meehan_all_fil_maf0.00_snps"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/stats_windows.py "$WDIR"/"$NAME".vcf 5000 2500 "$WDIR"/pop_structure.json "$WDIR"/outgroup_structure.json "$WDIR"/chr_list.txt
```

## Patterson's D (and associated)

Use biallelic SNP set to calculate D statistics for individuals and in a sliding window across genomes, using _E. grandis_ as the outgroup.

```bash
# Performed in UFRC queue system. See dsuite.job for more details.
# Resources:
module load gcc/9.3.0
module load dsuite/0.4-r49

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
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

Very few windows overlapping between outliers in the different statistics. I think this may be because the introgression is not fixed and therefore won't be detected by the more conservative way that I set cutoffs.

## Calculate genome scan statistics... again.

Instead, calculated pi for all pairwise comparisons between _individual samples_ of introgressants against reference populations.
```bash
# Performed in UFRC queue system. See pi_windows.job for more details.
# Resources: 200 Mb, 3 hrs

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pi_samplewise"

conda activate "$ENV_DIR"/euc_hyb_reseq

cd "$OUTDIR"

python "$SCRIPT_DIR"/pi_pairwise.py "$WDIR"/no_outgroup/meehan_all_fil_maf0.00_snps_biallelic_formatted.vcf 100000 20000 "$WDIR"/pop_structure.json "$WDIR"/outgroup_structure.json "$WDIR"/chr_list.txt
```

Calculate Dxy for Meehan Range _E. globulus_ versus _E. cordata_.

```bash
# Performed in UFRC queue system. See dxy_cordata.job for more details.
# Resources: 200 Mb, 6 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/dxy_cordata.py "$WDIR"/no_outgroup/meehan_all_fil_maf0.00_snps_biallelic_formatted.vcf 100000 20000 "$WDIR"/pop_structure.json "$WDIR"/outgroup_structure.json "$WDIR"/chr_list.txt
```

## Analysis of scan results... again.

Processed pi scan for each sample.
```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

ls *_scan_pi.tab | while read FILE
do
    Rscript "$SCRIPT_DIR"/process_pi_scans.r "$FILE"
done 2>&1 >> process_pi.out
```

Processed dxy scan against _E. cordata_ on local computer.
```R
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/windows_output/no_outgroup"
setwd(working_dir)

dxy_table_name <- "meehan_all_fil_maf0.00_snps_biallelic_formatted.vcf_dxy_cord.tab"
dxy_table <- read.table(dxy_table_name, header = TRUE, sep = "\t")

summary(dxy_table$Dxy) # Mean = 0.1953, Median = 0.1927
hist(dxy_table$Dxy) # Too tall to be normally distributed, but at least has the right shape

dxy_mean <- mean(dxy_table$Dxy, na.rm = TRUE)
dxy_sd <- sd(dxy_table$Dxy, na.rm = TRUE) # sd = 0.04419

flagged_dxy <- which(dxy_table$Dxy < dxy_mean - 2*dxy_sd)
flagged_dxy_coords <- dxy_table[flagged_dxy ,c("chr", "start", "end")] # 523 outlier windows

write.table(flagged_dxy_coords, "dxy_cord_flagged.tab", quote = FALSE, row.names = FALSE, col.names = TRUE)
```

Took a look at overlaps in windows of a few samples by hand before automating discovery.
```R
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/windows_output/no_outgroup"
setwd(working_dir)

WA01_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pi_samplewise/WA01_scan_pi.tab.flagged_pi.tab"
WA01_pi <- read.table(WA01_name, header = TRUE) # 799 outlier windows
intersect(flagged_dxy_coords, WA01_pi) # none

WA03_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pi_samplewise/WA03_scan_pi.tab.flagged_pi.tab" 
WA03_pi <- read.table(WA03_name, header = TRUE) # 796 outlier windows
intersect(flagged_dxy_coords, WA01_pi)

# ...etc. for each sample
```

No overlaps with windows identified as low distance to _E. cordata_.
NOTE: DID THIS WRONG, NEED TO TRY AGAIN WITH UPDATED CODE FOR FINDING COMMON WINDOWS.

## Genome scan with smaller windows
### 5k windows

Calculate all stats in 5kbp windows.

```bash
# Performed in UFRC queue system. See .job for more details.
# Resources:

```

```bash
module load R/4.2

SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/no_outgroup/5k_windows"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

# process stat windows to get outliers
cd "$WDIR"
Rscript "$SCRIPT_DIR"/reformat_fst.r all_maf0.00_no_outgr_5.windowed.weir.fst

Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_5k.tab Pi 2 above # mean: 16.09, sd = 15.24
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_5k.tab Deta 2 below
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_5k.tab Dxy 2 above
Rscript "$SCRIPT_DIR"/process_stat_scans.r all_maf0.00_no_outgr_5.windowed.weir.fst mean_FST 2 above
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_intr_cord_dxy_5k.tab Dxy 2 below

# shorten names for next part
rename 5k.tab 5k meehan_glob_pi_deta_dxy_5k.tab*
rename meehan_glob_pi_deta_dxy_5k 5k meehan_glob_pi_deta_dxy_5k*
rename .flagged_ _ 5k.flagged_*
mv all_maf0.00_no_outgr_5.windowed.weir.fst.flagged_mean_FST.tab 5k_FST.tab
mv meehan_intr_cord_dxy_5k.tab.flagged_Dxy.tab 5k_cord_Dxy.tab

# process outlier pi windows for each individual E. glob introgressant sample using pi mean and sd from all E. glob samples
cd "$WDIR"/pi_samplewise
while read NAME
do
    Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_pairwise_5k_"$NAME".tab Pi 2 above 16.09 15.24 2>&1 >> pi_samplewise_processing.out
done < "$LIST_DIR"/Eglobulus_MR.txt

# shorten names for next part
rename meehan_glob_pi_pairwise_5k_ 5k_Pi_ meehan_glob_pi_pairwise_5k_*.flagged_Pi.tab
rename .tab.flagged_Pi.tab .tab *

# Get overlaps between different outlier windows
ln -s ../5k_Deta.tab
ln -s ../5k_Dxy.tab
ln -s ../5k_FST.tab
ln -s ../5k_cord_Dxy.tab
printf "5k_Dxy.tab\n5k_FST.tab\n5k_cord_Dxy.tab\n" > overlap_files_1.txt
ls 5k_Pi_W*.tab > overlap_files_2.txt

Rscript "$SCRIPT_DIR"/compare_windows.r overlap_files_1.txt overlap_files_2.txt

# Then moved union files into new subdirectories for organization

# Get overlaps between common windows for pi and Dxy and Tajima's D
cd "$WDIR"/pi_samplewise/common_windows/Dxy_Pi
ln -s "$WDIR"/5k_Deta.tab
printf "5k_Deta.tab\n" > overlap_files_3.txt
ls common*Dxy*Pi*.tab > overlap_files_4.txt

Rscript "$SCRIPT_DIR"/compare_windows.r overlap_files_3.txt overlap_files_4.txt
rename common_5k_Deta.tab_common_5k_Dxy.tab_5k_Pi common_5k_Deta_Dxy_Pi *
# Moved Deta/Pi/Dxy overlap to its own subdirectory
```

### 10k windows

Identify outlier windows for each genome scan statistic and common outlier windows between stats.

```bash
module load R/4.2

SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/no_outgroup/10k_windows"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

# process stat windows to get outliers
cd "$WDIR"
Rscript "$SCRIPT_DIR"/reformat_fst.r all_maf0.00_no_outgr_10.windowed.weir.fst

Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_10k.tab Pi 2 above # mean: 28.71, sd: 26.2
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_10k.tab Deta 2 below
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_deta_dxy_10k.tab Dxy 2 above
Rscript "$SCRIPT_DIR"/process_stat_scans.r all_maf0.00_no_outgr_10.windowed.weir.fst mean_FST 2 above
Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_intr_cord_dxy_10k.tab Dxy 2 below

# shorten names for next part
rename 10k.tab 10k meehan_glob_pi_deta_dxy_10k.tab*
rename meehan_glob_pi_deta_dxy_10k 10k meehan_glob_pi_deta_dxy_10k*
rename .flagged_ _ 10k.flagged_*
mv all_maf0.00_no_outgr_10.windowed.weir.fst.flagged_mean_FST.tab 10k_FST.tab
mv meehan_intr_cord_dxy_10k.tab.flagged_Dxy.tab 10k_cord_Dxy.tab

# process outlier pi windows for each individual E. glob introgressant sample using pi mean and sd from all E. glob samples
cd "$WDIR"/pi_samplewise
while read NAME
do
    Rscript "$SCRIPT_DIR"/process_stat_scans.r meehan_glob_pi_pairwise_10k_"$NAME".tab Pi 2 above 28.71 26.2 2>&1 >> pi_samplewise_processing.out
done < "$LIST_DIR"/Eglobulus_MR.txt

# shorten names for next part
rename meehan_glob_pi_pairwise_10k_ 10k_Pi_ meehan_glob_pi_pairwise_10k_*.flagged_Pi.tab
rename .tab.flagged_Pi.tab .tab *

# Get overlaps between different outlier windows
ln -s ../10k_Deta.tab
ln -s ../10k_Dxy.tab
ln -s ../10k_FST.tab
ln -s ../10k_cord_Dxy.tab
printf "10k_Dxy.tab\n10k_FST.tab\n10k_cord_Dxy.tab\n" > overlap_files_1.txt
ls 10k_Pi_W*.tab > overlap_files_2.txt

Rscript "$SCRIPT_DIR"/compare_windows.r overlap_files_1.txt overlap_files_2.txt

# No common windows between outliers of pi and dxy, so stopped here.
```