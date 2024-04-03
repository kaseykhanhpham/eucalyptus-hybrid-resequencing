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

### Get genome-wide average of `pixy` results
Consolidated results.
```bash
# all but outgroup
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"

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
```

Analyzed results in `R` on local computer. Re-calculated pi and dxy genome-wide as suggested by `pixy` documentation.
```R
library(ggplot2)
library(dplyr)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy"
setwd(working_dir)

# import pixy files
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", na.strings="NA")
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", na.strings="NA")
all_fst <- read.table("all_fst.txt", header = TRUE, sep = "\t", na.strings="NA")

# filter by number of sites used to calculate in window
all_pi_fil <- all_pi[which(all_pi$no_sites > 40),]
all_dxy_fil <- all_dxy[which(all_dxy$no_sites > 40),]
all_fst_fil <- all_fst[which(all_fst$no_snps > 15),]

# PI
allpops_pi <- sum(all_pi_fil$count_diffs)/sum(all_pi_fil$count_comparisons) # 0.02827
glob_pi <- sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_comparisons"]) # 0.02790
ref_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02802
cord_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "cord"), "count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "cord"), "count_comparisons)"])

# DXY
allpops_dxy <- sum(all_dxy_fil$count_diffs)/sum(all_dxy_fil$count_comparisons) # 0.03554
ref_mr_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.02794
ref_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.04022
mr_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"])
```

Genome-wide results for pi:
| Statistic | Pop 1      | Value    |
| --------- | ---------- | -------- |
| pi        | all        | 0.02827  |
| pi        | glob (all) | 0.02790  |
| pi        | glob_MR    | 0.02788  |
| pi        | glob_pure  | 0.02802  |
| pi        | cord_MR    | 0.03014  |

Genome-wide results for dxy:
| Statistic | Pop 1     | Pop 2     | Value    |
| --------- | --------- | --------- | -------- |
| dxy       | all       | all       | 0.03554  |
| dxy       | glob_MR   | glob_pure | 0.02794  |
| dxy       | glob_MR   | cord_MR   | 0.04022  |
| dxy       | glob_pure | cord_MR   | 0.04044  |

### Graph `pixy` results
```R
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", as.is = TRUE)
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", as.is = TRUE)

# PI
## Meehan Range E. globulus
summary(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.glob MR", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.glob MR", xlab = "Log(Pi)")
abline(a = mean(log(all_pi[which(all_pi$pop == "glob_MR" & all_pi$no_sites >= 40), "avg_pi"])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00828 0.01734 0.02894 0.03750 0.38726

## Reference E. globulus
summary(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.glob ref", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "glob_pure" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.glob ref", xlab = "Log(Pi)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008038 0.017163 0.028810 0.037109 0.400848

## E. cordata
summary(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"])
hist(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"], breaks = 50, main = "Pi Distr E.cord", xlab = "Pi")
hist(log(all_pi[which(all_pi$pop == "cord_MR" & all_pi$no_sites >= 40), "avg_pi"]), breaks = 50, main = "Log(Pi) Distr E.cord", xlab = "Log(Pi)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008988 0.018740 0.031670 0.041366 0.392020

# DXY
## Meehan Range E. globulus & E. cordata
summary(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob MR - E.cord", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob MR - E.cord", xlab = "Log(Dxy)")
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01625 0.02917 0.04232 0.05646 0.38499

## Meehan Range E. globulus & Reference E. globulus
summary(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob MR - E.glob ref", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob MR - E.glob ref", xlab = "Log(Dxy)")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.008429 0.017479 0.028887 0.037468 0.386034

## Reference E. globulus & E. cordata
summary(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"])
hist(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"], breaks = 50, main = "Dxy Distr E.glob ref - E.cord", xlab = "Dxy")
hist(log(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR" & all_pi$no_sites >= 40), "avg_dxy"]), breaks = 50, main = "Log(Dxy) Distr E.glob ref - E.cord", xlab = "Log(Dxy)")
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.01607 0.02884 0.04192 0.05590 0.38714
```

## Patterson's D and associated statistics

Filtered SNP set to biallelic SNPs only using `vcftools`.

```bash
# Run on UFRC queue system; see get_biallelic.job for more details.
# Resources used: 10 Mb, 2 min

module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > "$VCFDIR"/all_fil_biallelic.vcf.gz # unzipped this file later for another analysis, so name will be different in files.
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

## Tajima's D and fixed differences
Calculated Tajima's D in 50,000 bp sliding windows using `vcftools`.

```bash
# Run on UFRC queue system; see tajd_windows.job for more details.
# Resources used: 5 Mb, 3 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Ecordata.txt --TajimaD 50000 --out cord_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_MR.txt --TajimaD 50000 --out glob_mr_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_ref.txt --TajimaD 50000 --out glob_ref_"$NAME"
done
```

Merged chromosome files.
```bash
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

touch cord_all.Tajima.D 
touch glob_mr_all.Tajima.D
touch glob_ref_all.Tajima.D

for NAME in "${VCFLIST[@]}"
do
    cat cord_"$NAME".Tajima.D >> cord_all.Tajima.D
    cat glob_mr_"$NAME".Tajima.D >> glob_mr_all.Tajima.D
    cat glob_ref_"$NAME".Tajima.D >> glob_ref_all.Tajima.D
done

rm *chr* # clean up chr-specific files
```

Calculated pi for individual sites for both reference populations using `vcftools`. I was looking for individual sites that have pi = 0 within species but pi != 0 between species as a signal of fixed differences. This excludes multi-allelic sites where one species has an allele fixed and the other has multiple other alleles unfixed, but given that biallelic SNPs make up ~1 million of the ~1.2 million SNPs total, I don't expect that particular scenario to seriously bias my results.

```bash
# Run on UFRC queue system; see pi_sites.job for more details.
# Resources used: 

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"

vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep Ecordata.txt --site-pi --out cord
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep Eglobulus_ref.txt --site-pi --out globref
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep ref_samples.txt --site-pi --out refs
```

Identified sites where pi = 0 within species and pi != 0 between species using `R`.
```R
# Import files
glob_tab <- read.table("globref.sites.pi", header = TRUE, sep = "\t")
cord_tab <- read.table("cord.sites.pi", header = TRUE, sep = "\t")
refs_tab <- read.table("refs.sites.pi", header = TRUE, sep = "\t")

# Create storage variables for loops
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "ChrUn")
chr_full <- c()
sites_full <- c()
for(chr in chr_list){
    # subset tables to just the current chromosome
    glob_fil <- glob_tab[which(glob_tab$CHROM == chr),]
    cord_fil <- cord_tab[which(cord_tab$CHROM == chr),]
    refs_fil <- refs_tab[which(refs_tab$CHROM == chr),]
    # Get full list of sites on the chromosome that appear in any of the tables
    chr_sites <- order(unique(c(glob_fil$POS, cord_fil$POS, refs_fil$POS)))
    sites_include <- c()
    for(site in chr_sites){
        # Check for presence of site in each table and if present, the value
        # looking for 0 in glob and cord only and != 0 in both table
        if(site %in% glob_fil$POS){
            glob_status <- glob_fil[which(glob_fil$POS == site), "PI"] == 0
        } else {
            glob_status <- FALSE
        }
        if(site %in% cord_fil$POS){
            cord_status <- cord_fil[which(cord_fil$POS == site), "PI"] == 0
        } else {
            cord_status <- FALSE
        }
        if(site %in% refs_fil$POS){
            refs_status <- refs_fil[which(refs_fil$POS == site), "PI"] != 0
        } else {
            refs_status <- FALSE
        }
        # add site to final list if 0 in glob and cord tables and != 0 in table with both
        if(glob_status & cord_status & refs_status){
            sites_include <- c(sites_include, site)
        }
    }
    # after looping through all sites on a chromosome, append the list to the genome-wide list
    chr_full <- c(chr_full, rep(chr, length(sites_include)))
    sites_full <- c(sites_full, sites_include)
}
# Write final list of fixed differences to output file
fdiffs_tab <- data.frame(chr = chr_full, pos = sites_full)
write.table(fdiffs_tab, "fixed_diffs_globref_cord.tab", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

See `chr_plots` README for synthesis of the results and plotting notes.