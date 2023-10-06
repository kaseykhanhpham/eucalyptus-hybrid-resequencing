# Sliding Window Genome Statistics

## Nucleotide Diversity and Divergence

Used [`pixy`](https://pixy.readthedocs.io/en/latest) to calculate pi (nucleotide diversity), dxy (absolute divergence), and FST (relative divergence) across the genome. First, removed the outgroup _E. grandis_ sample from VCFs and the outlier sample WF03/1051 for checking.

```bash
# Done on UFRC queue system; see remove_outgr_outl.job for more details.
# Resources used: ??? slurm fritzed out on me.

module load bcftools/1.15

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
OUTGR_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outgr"
OUTL_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outl"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    bcftools view --samples ^SRR10339635 --output-type z --output "$OUTGR_DIR"/"$NAME"_fil_no_outgr.vcf.gz --threads 12 "$VCF_DIR"/"$NAME"_fil.vcf.gz
    bcftools view --samples ^WF03,SRR10339635 --output-type z --output "$OUTL_DIR"/"$NAME"_fil_no_outl.vcf.gz --threads 12 "$VCF_DIR"/"$NAME"_fil.vcf.gz
done
```

Moved to `/orange` and indexed with `tabix`.

```bash
module load bcftools/1.15
OUTGR_DIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/06.snp_filtering/no_outgroup"
OUTL_DIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/06.snp_filtering/no_outlier"

cd "$OUTGR_DIR"
ls *.vcf.gz | while read FILE; do tabix "$FILE"; done

cd "$OUTL_DIR"
ls *.vcf.gz | while read FILE; do tabix "$FILE"; done
```

Created symlinks of all to `/blue` before proceeding.

Ran `pixy` for all samples except the _E. grandis_ outgroup.

```bash
# Done on UFRC queue system; see pixy.job for more details.
# Resources used: 1.3 Gb, 7 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outgr"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

conda activate "$ENV_DIR"/pixy

for NAME in "${VCFLIST[@]}"
do
    pixy --stats pi fst dxy --vcf "$VCF_DIR"/"$NAME"_fil_no_outgr.vcf.gz --populations "$LIST_DIR"/populations.txt --window_size 5000 --n_cores 12 --output_prefix "$NAME" --fst_type wc
done
```

| Statistic | Pop 1     | Pop 2     | Value    |
| --------- | --------- | --------- | -------- |
| pi        | glob_MR   | N/A       | 0.02842  |
| pi        | glob_pure | N/A       | 0.02856  |
| pi        | cord_MR   | N/A       | 0.03067  |
| dxy       | glob_MR   | glob_pure | 0.02847  |
| dxy       | glob_MR   | cord_MR   | 0.04104  |
| dxy       | glob_pure | cord_MR   | 0.04081  |

Ran `pixy` for all samples except WF03/1051 and the `E. grandis` outgroup.

```bash
# Done on UFRC queue system; see pixy_outlier.job for more details.
# Resources used: 1.3 Gb, 7 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outl"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/outlier_check"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

conda activate "$ENV_DIR"/pixy

for NAME in "${VCFLIST[@]}"
do
    pixy --stats pi fst dxy --vcf "$VCF_DIR"/"$NAME"_fil_no_outl.vcf.gz --populations "$LIST_DIR"/populations_no_outl.txt --window_size 5000 --n_cores 12 --output_prefix "$NAME" --fst_type wc
done
```

### Compare `pixy` results with and without outlier WF03/1051
Consolidated results:

```bash
# all but outgroup
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
OUTL_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/outl_check"
declare -a VCFLIST=(chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

cd "$WDIR"
cat chr01_pi.txt > all_pi.txt
cat chr01_dxy.txt > all_dxy.txt
cat chr01_fst.txt > all_fst.txt

for NAME in "${VCFLIST[@]}"
do
    tail -n +1 "$NAME"_pi.txt >> all_pi.txt
    tail -n +1 "$NAME"_dxy.txt >> all_dxy.txt
    tail -n +1 "$NAME"_fst.txt >> all_fst.txt
done

cd "$OUTL_DIR"
cat chr01_pi.txt > all_outl_pi.txt
cat chr01_dxy.txt > all_outl_dxy.txt
cat chr01_fst.txt > all_outl_fst.txt

for NAME in "${VCFLIST[@]}"
do
    tail -n +1 "$NAME"_pi.txt >> all_outl_pi.txt
    tail -n +1 "$NAME"_dxy.txt >> all_outl_dxy.txt
    tail -n +1 "$NAME"_fst.txt >> all_outl_fst.txt
done
```

Analyzed results in `R` on local computer.
```R
library(ggplot2)
library(dplyr)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy"
setwd(working_dir)

# import pixy files
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", na.strings="NA")
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", na.strings="NA")
all_fst <- read.table("all_fst.txt", header = TRUE, sep = "\t", na.strings="NA")

outl_pi <- read.table("all_outl_pi.txt", header = TRUE, sep = "\t", na.strings="NA")
outl_dxy <- read.table("all_outl_dxy.txt", header = TRUE, sep = "\t", na.strings="NA")
outl_fst <- read.table("all_outl_fst.txt", header = TRUE, sep = "\t", na.strings="NA")

# filter by number of sites used to calculate in window
all_pi_fil <- all_pi[which(all_pi$no_sites > 40),]
outl_pi_fil <- outl_pi[which(outl_pi$no_sites > 40),]

all_dxy_fil <- all_dxy[which(all_dxy$no_sites > 40),]
outl_dxy_fil <- outl_dxy[which(outl_dxy$no_sites > 40),]

all_fst_fil <- all_fst[which(all_fst$no_snps > 15),]
outl_fst_fil <- outl_fst[which(outl_fst$no_snps > 15),]

# compare genome-wide values of pi and dxy for datasets with and without the outlier, WF03/1051

# PI
allpops_pi <- sum(all_pi_fil$count_diffs)/sum(all_pi_fil$count_comparisons) # 0.02827
allpops_pi_outl <- sum(outl_pi_fil$count_diffs)/sum(outl_pi_fil$count_comparisons) # 0.02843

glob_pi <- sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_comparisons"]) # 0.02790
glob_pi_outl <- sum(outl_pi_fil[which(outl_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_diffs"])/sum(outl_pi_fil[which(outl_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_comparisons"]) # 0.02809

ref_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02802
ref_pi_outl <- sum(outl_pi_fil[which(outl_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(outl_pi_fil[which(outl_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02938

# DXY
allpops_dxy <- sum(all_dxy_fil$count_diffs)/sum(all_dxy_fil$count_comparisons) # 0.03554
allpops_dxy_outl <- sum(outl_dxy_fil$count_diffs)/sum(outl_dxy_fil$count_comparisons) # 0.03610

ref_mr_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.02794
ref_mr_dxy_outl <- sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_MR" & outl_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_MR" & outl_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.02857

ref_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.04022
ref_cord_dxy_outl <- sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_pure" & outl_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_pure" & outl_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.04089

mr_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"])
```

Whole-genome pi comparisons:
| group      | pi with WF03/1051 | pi without WF03/1051 |
| ---------- | ----------------- | -------------------- |
| all        | 0.02827           | 0.02843              |
| _E. glob_  | 0.02790           | 0.02809              |
| _glob_ ref | 0.02802           | 0.02938              |
| _glob_ MR  | 0.02788           | N/A                  |
| _cord_     | 0.03014           | N/A                  |

Whole-genome dxy comparisons:
| group 1   | group 2    | dxy with WF03/1051 | dxy without WF03/1051 |
| --------- | ---------- | ------------------ | --------------------- |
| all       | all        | 0.03554            | 0.03610               |
| _glob_ MR | _glob_ ref | 0.02794            | 0.02857               |
| _cord_    | _glob_ ref | 0.04022            | 0.04089               |
| _glob_ MR | _cord_     | 0.04044            | N/A                   |

Inclusion of WF03/1051 changed the estimated genome-wide pi for the reference group alone the most, as expected. Values of dXY remained mostly not too strongly different. The numbers don't seem so extremely different that I would be warranted in excluding the sample...

### Characterize `pixy` results
```R
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", as.is = TRUE)
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", as.is = TRUE)

summary(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.glob MR", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.glob MR", xlab = "Log(Pi)")
abline(a = mean(log(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00828 0.01734 0.02894 0.03750 0.38726

summary(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.glob ref", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.glob ref", xlab = "Log(Pi)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008038 0.017163 0.028810 0.037109 0.400848
summary(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.cord", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.cord", xlab = "Log(Pi)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008988 0.018740 0.031670 0.041366 0.392020

summary(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob MR - E.cord", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob MR - E.cord", xlab = "Log(Dxy)")
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01625 0.02917 0.04232 0.05646 0.38499
summary(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob MR - E.glob ref", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob MR - E.glob ref", xlab = "Log(Dxy)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008429 0.017479 0.028887 0.037468 0.386034
summary(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob ref - E.cord", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob ref - E.cord", xlab = "Log(Dxy)")
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01607 0.02884 0.04192 0.05590 0.38714
```

### Identify outlier windows
Split pixy output files by taxon/comparison.

```R
options(scipen=999) 

all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", as.is = TRUE)
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", as.is = TRUE)

## recalculate pi for all E. globulus samples
get_pi <- function(df){
    # recalculate pi
    tab_subset <- df[which(df$pop %in% c("glob_MR", "glob_pure")),]
    new_count_diffs <- sum(as.numeric(tab_subset$count_diffs), na.rm = TRUE)
    new_count_comparisons <- sum(as.numeric(tab_subset$count_comparisons), na.rm = TRUE)
    new_avg_pi <- new_count_diffs/new_count_comparisons

    # re-calibrate window positions to be in standard genome coords (0-index)
    new_window_pos_1 <- as.numeric(tab_subset[1, "window_pos_1"]) - 1
    new_window_pos_2 <- as.numeric(tab_subset[1, "window_pos_2"]) - 1
    
    # sum numbers for numeric fields
    new_count_missing <- sum(as.numeric(tab_subset$count_missing), na.rm = TRUE)
    new_no_sites <- sum(as.numeric(tab_subset$no_sites))

    # export processed row (convert all vars back to characters for rbind)
    new_row <- c("glob_all", tab_subset[1, "chromosome"], as.character(new_window_pos_1), as.character(new_window_pos_2), as.character(new_avg_pi), as.character(new_no_sites), as.character(new_count_diffs), as.character(new_count_comparisons), as.character(new_count_missing))
    return(new_row)
}

new_glob_rows <- sapply(c(1:length(which(all_pi$pop == "glob_MR"))), function(x) get_pi(all_pi[(x*3-2):(x*3),]))
new_glob_tab <- t(new_glob_rows)
colnames(new_glob_tab) <- c("pop", "chromosome", "window_pos_1", "window_pos_2", "avg_pi", "no_sites", "count_diffs", "count_comparisons", "count_missing")
write.table(new_glob_tab, "glob_pi_formatted.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

## process dxy for E. cordata - E. globulus MR comparisons
subset_dxy <- all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure"),]

# adjust window positions
subset_dxy$window_pos_1 <- as.numeric(subset_dxy$window_pos_1) - 1
subset_dxy$window_pos_2 <- as.numeric(subset_dxy$window_pos_2) - 1

write.table(subset_dxy, "glob_MR_cord_MR_dxy.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
```

Retrieved outlier windows using custom `R` script.

```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_pi_formatted.txt sd 2,3,4,6,5 2 above 40 glob_pi_outliers_sd2.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_pi_formatted.txt percent 2,3,4,6,5 0.05 above 40 glob_pi_outliers_p05.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_pi_formatted.txt percent 2,3,4,6,5 0.10 above 40 glob_pi_outliers_p10.txt

Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_MR_cord_MR_dxy.txt sd 3,4,5,7,6 2 below 40 glob_MR_cord_MR_dxy_outliers_sd2.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_MR_cord_MR_dxy.txt percent 3,4,5,7,6 0.05 below 40 glob_MR_cord_MR_dxy_outliers_p95.txt

Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_MR_cord_MR_dxy.txt percent 3,4,5,7,6 0.10 below 40 glob_MR_cord_MR_dxy_outliers_p90.txt
```

| Stat  | Taxon        | Threshold    | Number of Windows |
| ----- | ------------ | ------------ | ----------------- |
| pi    | glob         | > +2 sd      | 3049              |
| pi    | glob         | top 5%       | 2743              |
| pi    | glob         | top 10%      | 4738              |
| dxy   | glob MR/cord | < -2 sd      | 0                 |
| dxy   | glob MR/cord | bottom 5%    | 2743              |
| dxy   | glob MR/cord | bottom 10%   | 4738              |

Converted `pixy` output into BED file of windows.

```R
wdir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
setwd(wdir)
options(scipen=999)

all_dxy <- read.table("all_dxy.txt", header = TRUE)
all_pi <- read.table("all_pi.txt", header = TRUE)

dxy_bed <- all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
dxy_bed$window_pos_1 <- as.numeric(dxy_bed$window_pos_1) - 1
write.table(dxy_bed, "all_dxy.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

pi_bed <- all_pi[which(all_pi$pop == "glob_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
pi_bed$window_pos_1 <- as.numeric(pi_bed$window_pos_1) - 1
write.table(pi_bed, "all_pi.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
```

## Patterson's D and associated statistics

Filtered SNP set to biallelic SNPs only using `vcftools`.

```bash
# Run on UFRC queue system; see get_biallelic.job for more details.
# Resources used: 10 Mb, 2 min

module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > "$VCFDIR"/all_fil_biallelic.vcf.gz
```

Calculated fD and fDM in sliding windows of 40 viable SNPs at a time using [`Dsuite`](https://github.com/millanek/Dsuite).

```bash
# Run on UFRC queue system; see dsuite.job for more details.
# Resources used: 4 Mb, 8 min

module load gcc/12.2.0

DSUITE_DIR="/blue/soltis/kasey.pham/bin/Dsuite/Build"
# My install is compiled on GCC 12.2.0, versus the HPC install on 9.3.0
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"

"$DSUITE_DIR"/Dsuite Dinvestigate -w 40,20 -g "$VCF_DIR"/all_fil_biallelic.vcf.gz SETS.txt test_trios.txt
mv glob_pure_glob_MR_cord_MR_localFstats__40_20.txt localFstats_40_20.txt
```

### Get outlier windows
Retrieved outlier windows using custom `R` script.

```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

# get highest 5% of windows for each sliding window statistic
Rscript "$SCRIPT_DIR"/process_stat_windows.r localFstats_40_20.txt percent 1,2,3,-1,5 0.5 above 0 fD_40_20_outliers_p05.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r localFstats_40_20.txt percent 1,2,3,-1,6 0.5 above 0 fDm_40_20_outliers_p05.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r localFstats_40_20.txt percent 1,2,3,-1,5 0.5 above 0 df_40_20_outliers_p05.txt
```


## Tajima's D and Transition/Transversion
Calculated Tajima's D and Ts/Tv rate in 5000bp sliding windows using `vcftools`.

```bash
# Run on UFRC queue system; see tajd_windows.job for more details.
# Resources used: 5 Mb, 3 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Ecordata.txt --TajimaD 5000 --out cord_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_MR.txt --TajimaD 5000 --out glob_mr_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_ref.txt --TajimaD 5000 --out glob_ref_"$NAME"
done
```

```bash
# Run on UFRC queue system; see tstv_windows.job for more details.
# Resources used: 5 Mb, 2 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Ecordata.txt --TsTv 5000 --out cord_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_MR.txt --TsTv 5000 --out glob_mr_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_ref.txt --TsTv 5000 --out glob_ref_"$NAME"
done
```

Merged chromosome files.

```bash
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

touch cord_all.Tajima.D 
touch cord_all.TsTv
touch glob_mr_all.Tajima.D
touch glob_mr_all.TsTv
touch glob_ref_all.Tajima.D
touch glob_ref_all.TsTv

for NAME in "${VCFLIST[@]}"
do
    cat cord_"$NAME".Tajima.D >> cord_all.Tajima.D
    cat cord_"$NAME".TsTv >> cord_all.TsTv
    cat glob_mr_"$NAME".Tajima.D >> glob_mr_all.Tajima.D
    cat glob_mr_"$NAME".TsTv >> glob_mr_all.TsTv
    cat glob_ref_"$NAME".Tajima.D >> glob_ref_all.Tajima.D
    cat glob_ref_"$NAME".TsTv >> glob_ref_all.TsTv
done

rm *chr* # clean up chr-specific files
```

### Get outlier windows

Retrieved outlier windows using custom `R` script.

```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_mr_all.Tajima.D percent 1,2,-1,3,4 0.05 below 40 glob_mr_all_tajd_outliers_low.txt
Rscript "$SCRIPT_DIR"/process_stat_windows.r glob_mr_all.Tajima.D percent 1,2,-1,3,4 0.05 above 40 glob_mr_all_tajd_outliers_high.txt
```

## Compare outlier windows

```R
library(dplyr)

# import tables of outliers for each statistic
fd_outl <- read.table("./dsuite/fD_40_20_outliers_p05.txt", header = TRUE)
fdm_outl <- read.table("./dsuite/fDm_40_20_outliers_p05.txt", header = TRUE)
df_outl <- read.table("./dsuite/df_40_20_outliers_p05.txt", header = TRUE)

dxy_outl_p95 <- read.table("./pixy/glob_MR_cord_MR_dxy_outliers_p95.txt", header = TRUE)
dxy_outl_p90 <- read.table("./pixy/glob_MR_cord_MR_dxy_outliers_p90.txt", header = TRUE)
dxy_outl_p85 <- read.table("./pixy/glob_MR_cord_MR_dxy_outliers_p85.txt", header = TRUE)
dxy_outl_sd2 <- read.table("./pixy/glob_MR_cord_MR_dxy_outliers_sd2.txt", header = TRUE)
pi_outl_p05 <- read.table("./pixy/glob_pi_outliers_p05.txt", header = TRUE)
pi_outl_p10 <- read.table("./pixy/glob_pi_outliers_p10.txt", header = TRUE)
pi_outl_p15 <- read.table("./pixy/glob_pi_outliers_p15.txt", header = TRUE)
pi_outl_sd2 <- read.table("./pixy/glob_pi_outliers_sd2.txt", header = TRUE)

tajd_outl <- read.table("./vcftools/glob_mr_all_tajd_outliers.txt", header = TRUE)

# compare windows for overlaps
dxy05_pi05 <- inner_join(dxy_outl_p95, pi_outl_p05) # 3 common
dxy10_pi10 <- inner_join(dxy_outl_p90, pi_outl_p10) # 11 common
dxy15_pi15 <- inner_join(dxy_outl_p85, pi_outl_p15) # 22 common

write.table(dxy05_pi05, "pi_dxy_outl_p05.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(dxy10_pi10, "pi_dxy_outl_p10.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(dxy15_pi15, "pi_dxy_outl_p15.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
```

Removed headers from Dxy/Pi overlap files to make them BED format

```bash
PIXY_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
declare -a PIXY_LIST=(pi_dxy_outl_p05 pi_dxy_outl_p10 pi_dxy_outl_p15)
DSUITE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"
declare -a DSUITE_LIST=(df_40_20_outliers_p05 fD_40_20_outliers_p05 fDm_40_20_outliers_p05)

for NAME in "${PIXY_LIST[@]}"
do
    tail -n +2 "$PIXY_DIR"/"$NAME".txt | awk '{gsub(" ","\t"); print}' > "$PIXY_DIR"/"$NAME".bed
done

for NAME in "${DSUITE_LIST[@]}"
do
    tail -n +2 "$DSUITE_DIR"/"$NAME".txt | awk '{gsub(" ","\t"); print}' > "$DSUITE_DIR"/"$NAME".bed
done
```

Get overlap in fdM windows and pi/dxy windows using `BEDTools`.

```bash
PIXY_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
DSUITE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"

module load bedtools/2.30.0

bedtools intersect -a "$PIXY_DIR"/pi_dxy_outl_p05.bed -b "$DSUITE_DIR"/df_40_20_outliers_p05.bed -wa > pi_dxy_outl_p05_df.bed
```

Fstats of overlapping windows as compared to Fstat distribution:

**5%**
| Chr   | Start    | End       | fD       | fDM      | df       |
| ----- | -------- | --------- | -------- | -------- | -------- |
| Chr08 | 67870000 | 67874999  | 0.127294 | 0.077470 | 0.007757 |

**10%**
| Chr   | Start    | End      | fD        | fDM       | df        |
| ----- | ---------| ---------| --------- | --------- | --------- |
| Chr02 | 8340000  | 8344999  | -0.006772 | -0.005344 | -0.008166 |
| Chr02 | 14115000 | 14119999 | 0.371875  | 0.278639  | 0.337718  |
| Chr03 | 17195000 | 17199999 | -0.053282 | -0.042062 | -0.029466 |
| Chr03 | 30380000 | 30384999 | 0.003036  | 0.002502  | 0.015761  |
| Chr03 | 57880000 | 57884999 | 0.051417  | 0.038933  | 0.028908  |
| Chr03 | 61200000 | 61204999 | 0.027935  | 0.022333  | 0.003827  |
| Chr05 | 59150000 | 59154999 | -0.015417 | -0.010327 | -0.022913 |
| Chr06 | 16120000 | 16124999 | -0.008713 | -0.007472 | -0.005123 |
| Chr06 | 23855000 | 23859999 | -0.280654 | -0.159607 | -0.051590 |
| Chr09 | 33970000 | 33974999 | -0.040472 | -0.026516 | -0.009094 |

**15%**
| Chr   | Start    | End      | fD         | fDM        | df        |
| ----- | ---------| ---------| ---------- | ---------- | --------- |
| Chr03 | 29230000 | 29234999 |  0.032828  | 0.025977   | 0.004450  |
| Chr03 | 61545000 | 61549999 | 0.042820   | 0.035960   | 0.005258  |
| Chr06 | 4010000  | 4014999  | 0.104500   | 0.075476   | 0.015816  |
| Chr06 | 17415000 | 17419999 |  -0.008713 | -0.007472  | -0.005123 |
| Chr06 | 23955000 | 23959999 | -0.280488  | -0.153720  | -0.041347 |
| Chr06 | 23975000 | 23979999 | -0.280488  | -0.153720  | -0.041347 |
| Chr06 | 24765000 | 24769999 | -0.028896  | -0.021621  | -0.022063 |
| Chr06 | 46335000 | 46339999 | 0.204911   | 0.171211   | 0.142243  |
| Chr08 | 6195000  | 6199999  | 0.039659   | 0.028836   | 0.019434  |
| Chr08 | 17045000 | 17049999 | -0.004604  | -0.002842  | 0.001133  |
| Chr09 | 36560000 | 36564999 |  -0.040472 | -0.026516  | -0.009094 |

## Calculate dXY for admixed individuals one-by-one

Created population files for each admixed individual.
```bash
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals"

while read NAME
do
    cat "$WDIR"/populations_cord.txt > "$WDIR"/population_files/populations_cord_"$NAME".txt
    printf '%s\tglob_MR\n' "$NAME" >> "$WDIR"/population_files/populations_cord_"$NAME".txt
done < "$WDIR"/Eglobulus_MR.txt
```

Ran `pixy` on all individual _E. globulus_ MR samples.

```bash
# Performed in UFRC queue system; see pixy_inds.job for more details.
# Resources used: 1.4 Gb, 2 hrs

module load conda
module load vcftools/0.1.16

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals"

conda activate "$ENV_DIR"/pixy

while read NAME
do
    pixy --stats dxy --vcf "$VCF_DIR"/all_fil.vcf.gz --populations "$LIST_DIR"/population_files/populations_cord_"$NAME".txt --window_size 5000 --n_cores 12 --output_prefix cord_"$NAME"
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Retrieved outlier windows for each sample.

```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals"

while read NAME
do
    Rscript "$SCRIPT_DIR"/process_stat_windows.r "$LIST_DIR"/cord_"$NAME"_dxy.txt percent 3,4,5,7,6 0.05 below 40 cord_"$NAME"_dxy_outl_p95.txt
done < "$LIST_DIR"/Eglobulus_MR.txt

while read NAME
do
    Rscript "$SCRIPT_DIR"/process_stat_windows.r "$LIST_DIR"/cord_"$NAME"_dxy.txt percent 3,4,5,7,6 0.10 below 40 cord_"$NAME"_dxy_outl_p90.txt
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Converted outlier windows to genome annotation BED files.

```bash
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals/outlier_files"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals"

while read NAME
do
    tail -n +2 cord_"$NAME"_dxy_outl_p90.txt | awk '{gsub(" ","\t"); print}' > "$WDIR"/bedfiles/cord_"$NAME"_dxy_outl_p90.bed
    tail -n +2 cord_"$NAME"_dxy_outl_p95.txt | awk '{gsub(" ","\t"); print}' > "$WDIR"/bedfiles/cord_"$NAME"_dxy_outl_p95.bed
done < "$LIST_DIR"/Eglobulus_MR.txt
```