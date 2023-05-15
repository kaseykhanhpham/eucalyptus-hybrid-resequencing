# Genome Scan Analysis

## Calculate statistics

### Nucleotide diversity, absolute divergence, selection
Calculated diagnostic statistics for introgression in sliding windows across genome using the python library `egglib`. Wrote a script to automate egglib calculations for sliding windows: `stats_windows.py`. It requires the following input files present in this directory:

`chr_list.txt`: list of chromosome/contig names in reference genome, one per line

**structure files**: The JSON files in this directory define the individuals and populations to be used for calculations of diagnostic stats. The "glob" files consider introgressed and reference _E. globulus_ to be the ingroup (stats calculated) and _E. cordata_ and _E. grandis_ to be the outgroup (stats not calculated). The "cord" files consider introgressed _E. globulus_ and _E. cordata_ to be the ingroup and reference _E. globulus_ and _E. grandis_ to be the outgroup.

Ran `stats_windows.py` with the stats arguments "Pi", "Dxy", "FST", "Deta" (see below).

```bash
# Done on UFRC queue system, see pi_windows.job, dxy_windows.job, dxy_cord_windows.job, fst_windows.job, deta_windows.job for more detail
# Code below is example from pi_windows.job
# Maximum resources used: 370 Mb, 35 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/stats_windows.py "$VCF_DIR"/meehan_all_fil_maf0.00_snps.vcf Pi 5000 2500 "$WDIR"/glob_structure.json "$WDIR"/glob_output.json "$WDIR"/chr_list.txt glob_pi_windows.tab
```

### Patterson's D (and associated)

First filtered VCF to only biallelic or monoallelic sites.

```bash
# Performed in UFRC queue system. See get_biallelic.job for more details.
# Resources: 4 Mb, 17 min
module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"

vcftools --vcf "$VCFDIR"/meehan_all_fil_maf0.00_snps.vcf --min-alleles 2 --max-alleles 2 --recode --stdout > meehan_all_fil_maf0.00_snps_biallelic.vcf
```

Used biallelic SNP set to calculate D statistics for individuals and in a sliding window across genomes, using _E. grandis_ as the outgroup.

```bash
# Performed in UFRC queue system. See dsuite.job for more details.
# Resources: 4 Mb, 17 min
DSUITE_DIR="/blue/soltis/kasey.pham/bin/Dsuite/Build"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"
NAME="meehan_all_fil_maf0.00_snps_biallelic"

"$DSUITE_DIR"/Dsuite Dinvestigate -n "$NAME" "$WDIR"/"$NAME".vcf "$WDIR"/SETS.txt test_trios.txt

# Moved Dsuite output to shorter name: glob-pure-cord_localFstats.txt
```

## Process outlier windows

Identified regions where statistic was +/- 2sd from the mean value across all windows using a custom `R` script, where stat, cutoff, and number of SNPs for inclusion can be specified.

```bash
module load R/4.2

SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

# Pi
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_pi_windows.tab Pi sd 2 above 25 
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_pi_windows.tab Pi percent 0.05 above 25 
# Dxy - E. globulus MR to E. globulus ref
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_dxy_windows.tab Dxy sd 2 above 25
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_dxy_windows.tab Dxy percent 0.05 above 25
# Dxy - E. globulus MR to E. cordata
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/cord_dxy_windows.tab Dxy sd 2 below 25
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/cord_dxy_windows.tab Dxy percent 0.05 below 25
# Tajima's D
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_tajd_windows.tab D sd 2 below 25
Rscript "$SCRIPT_DIR"/process_stat_windows.r "$WDIR"/glob_tajd_windows.tab D percent 0.05 below 25
```

Examined overlap between outlier windows in `R`:

```R
# done locally
library(dplyr)
# read in outlier tables
cord_dxy_sd <- read.table("cord_dxy_windows.tab.sd2.flagged.tab", header = TRUE, sep = " ") # 1179
glob_deta_sd <- read.table("glob_deta_windows.tab.sd2.flagged.tab", header = TRUE, sep = " ") # 1867
glob_dxy_sd <- read.table("glob_dxy_windows.tab.sd2.flagged.tab", header = TRUE, sep = " ") # 4809 
glob_pi_sd <- read.table("glob_pi_windows.tab.sd2.flagged.tab", header = TRUE, sep = " ") # 4718

cord_dxy_glob_pi_sd <- inner_join(cord_dxy_sd, glob_pi_sd) # 0
glob_dxy_glob_pi_sd <- inner_join(glob_dxy_sd, glob_pi_sd) # 4278
cord_dxy_glob_deta_sd <- inner_join(cord_dxy_sd, glob_deta_sd) # 270

cord_dxy_p <- read.table("cord_dxy_windows.tab.percent0.05.flagged.tab", header = TRUE, sep = " ") # 6417
glob_deta_p <- read.table("glob_deta_windows.tab.percent0.05.flagged.tab", header = TRUE, sep = " ") # 6417
glob_dxy_p <- read.table("glob_dxy_windows.tab.percent0.05.flagged.tab", header = TRUE, sep = " ") # 6418
glob_pi_p <- read.table("glob_pi_windows.tab.percent0.05.flagged.tab", header = TRUE, sep = " ") # 6418

cord_dxy_glob_pi_p <- inner_join(cord_dxy_p, glob_pi_p) # 0
glob_dxy_glob_pi_p <- inner_join(glob_dxy_p, glob_pi_p) # 5837
cord_dxy_glob_deta_p <- inner_join(cord_dxy_p, glob_deta_p) # 1898
```

## Calculate Pi, Dxy, FST with Pixy

Remove outgroup from VCF and zip the file.

```bash
module load vcftools
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"

vcftools --vcf "$VCF_DIR"/meehan_all_fil_maf0.00_snps.vcf --remove-indv SRR10339635 --recode --out "$WDIR"/meehan_all_fil_maf0.00_snps_noout.vcf
bgzip meehan_all_fil_maf0.00_snps_noout.vcf
tabix meehan_all_fil_maf0.00_snps_noout.vcf
```

Run pixy.

```bash

```