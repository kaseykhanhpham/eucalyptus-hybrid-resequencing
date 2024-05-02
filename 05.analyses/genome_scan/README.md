# Sliding Window Genome Statistics
## Nucleotide Diversity and Divergence
Used [`pixy`](https://pixy.readthedocs.io/en/latest) to calculate pi (nucleotide diversity), dxy (absolute divergence), and FST (relative divergence) across the genome. 

Created symlinks of all chromosome variant call VCFs and tabix indices to `/blue` for use by `pixy`.

Ran `pixy` for all samples.

```bash
# Done on UFRC queue system; see pixy.job for more details.
# Resources used: 1.3 Gb, 17 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

conda activate "$ENV_DIR"/pixy

for NAME in "${VCFLIST[@]}"
do
    pixy --stats pi fst dxy --vcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --populations "$LIST_DIR"/populations.txt --window_size 5000 --n_cores 12 --output_prefix "$NAME" --fst_type wc
done
```

### Get genome-wide average of `pixy` results
Consolidated results.
```bash
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
declare -a VCFLIST=(chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

cd "$WDIR"
cat chr01_pi.txt > all_pi.txt
cat chr01_dxy.txt > all_dxy.txt
cat chr01_fst.txt > all_fst.txt

for NAME in "${VCFLIST[@]}"
do
    tail -n +2 "$NAME"_pi.txt >> all_pi.txt
    tail -n +2 "$NAME"_dxy.txt >> all_dxy.txt
    tail -n +2 "$NAME"_fst.txt >> all_fst.txt
done

# clean-up
rm chr*.txt
```

Analyzed results in `R` on local computer. Re-calculated pi and dxy genome-wide as suggested by `pixy` documentation.
```R
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy"
setwd(working_dir)

# import pixy files
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", na.strings = "NA")
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", na.strings = "NA")
all_fst <- read.table("all_fst.txt", header = TRUE, sep = "\t", na.strings = "NA")

# unfiltered stats
# PI
sum(all_pi$count_diffs, na.rm = TRUE)/sum(all_pi$count_comparisons, na.rm = TRUE) # ALL: 0.02553995
sum(all_pi[which(all_pi$pop == "glob_MR"),"count_diffs"], na.rm = TRUE)/sum(all_pi[which(all_pi$pop == "glob_MR"),"count_comparisons"], na.rm = TRUE) # MR GLOB: 0.02523312
sum(all_pi[which(all_pi$pop == "glob_pure"),"count_diffs"], na.rm = TRUE)/sum(all_pi[which(all_pi$pop == "glob_pure"),"count_comparisons"], na.rm = TRUE) # REF GLOB: 0.02534997
sum(all_pi[which(all_pi$pop == "cord_MR"), "count_diffs"], na.rm = TRUE)/sum(all_pi[which(all_pi$pop == "cord_MR"), "count_comparisons"], na.rm = TRUE) # CORD: 0.02699768

# DXY
sum(all_dxy$count_diffs, na.rm = TRUE)/sum(all_dxy$count_comparisons, na.rm = TRUE) # 0.03227452
sum(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure"),"count_diffs"], na.rm = TRUE)/sum(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "glob_pure"),"count_comparisons"], na.rm = TRUE) # 0.0253001
sum(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR"),"count_diffs"], na.rm = TRUE)/sum(all_dxy[which(all_dxy$pop1 == "glob_pure" & all_dxy$pop2 == "cord_MR"),"count_comparisons"], na.rm = TRUE) # 0.03655834
sum(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR"),"count_diffs"], na.rm = TRUE)/sum(all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR"),"count_comparisons"], na.rm = TRUE) # 0.03677519

# filter by number of sites used to calculate in window
all_pi_fil <- all_pi[which(all_pi$no_sites > 40),]
all_dxy_fil <- all_dxy[which(all_dxy$no_sites > 40),]
all_fst_fil <- all_fst[which(all_fst$no_snps > 15),]

# PI
allpops_pi <- sum(all_pi_fil$count_diffs)/sum(all_pi_fil$count_comparisons) # 0.02518261
mr_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "glob_MR"), "count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "glob_MR"),"count_comparisons"]) # 0.02487323
ref_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02498456
cord_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "cord_MR"), "count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "cord_MR"), "count_comparisons"]) # 0.02665737

# DXY
allpops_dxy <- sum(all_dxy_fil$count_diffs)/sum(all_dxy_fil$count_comparisons) # 0.2062443
ref_mr_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.01000927
mr_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.03641528
ref_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.03620112

# FST
allpops_fst <- mean(all_fst_fil$avg_wc_fst, na.rm = TRUE) # 0.2062146
ref_mr_fst <- mean(all_fst_fil[which(all_fst_fil$pop1 == "glob_MR" & all_fst_fil$pop2 == "glob_pure"),"avg_wc_fst"], na.rm = TRUE) # 0.009491992
mr_cord_fst <- mean(all_fst_fil[which(all_fst_fil$pop1 == "glob_MR" & all_fst_fil$pop2 == "cord_MR"),"avg_wc_fst"], na.rm = TRUE) # 0.3079715
ref_cord_fst <- mean(all_fst_fil[which(all_fst_fil$pop1 == "glob_pure" & all_fst_fil$pop2 == "cord_MR"),"avg_wc_fst"], na.rm = TRUE) # 0.3007242
```

Genome-wide results for pi:
| Statistic | Pop 1      | Value      |
| --------- | ---------- | ---------- |
| pi        | all        | 0.02518261 |
| pi        | glob_MR    | 0.02487323 |
| pi        | glob_pure  | 0.02498456 |
| pi        | cord_MR    | 0.02665737 |

Genome-wide results for divergence:
| Pop 1     | Pop 2     | dxy Value  | FST Value  |
| --------- | --------- | ---------- | ---------- |
| all       | all       | 0.0319177  | 0.2062443  |
| glob_MR   | glob_pure | 0.02494443 | 0.01000927 |
| glob_MR   | cord_MR   | 0.03641528 | 0.3079715  |
| glob_pure | cord_MR   | 0.03620112 | 0.3007242  |

### Graph `pixy` results
In the same `R` environment as above. Uses code from [the R Graphics Cookbook](https://r-graphics.org/recipe-distribution-multi-hist).

```R
library(ggplot2)
library(dplyr)

# PI
all_pi_fil$pop_clean <- recode_factor(all_pi_fil$pop, "cord_MR" = "E. cordata", "glob_MR" = "MR  E. globulus", "glob_pure" = "Pure  E. globulus")
pi_means <- data.frame(pop_clean = c("E. cordata", "MR  E. globulus", "Pure  E. globulus"), mean = c(cord_pi, mr_pi, ref_pi))
pi_plot <- ggplot(all_pi_fil, aes(x = avg_pi, fill = pop_clean)) + geom_histogram(color = "black") + 
           scale_fill_manual(name = "", values = c("#FFA716", "#07adcd", "gray")) +
           facet_grid(pop_clean ~ .) + theme_bw() + labs(x = "Average Pi per 5kb Window", y = "Count") +
           geom_vline(data=pi_means, aes(xintercept = mean), colour="black") + 
           theme(legend.position = "none", axis.title=element_text(size=14))

# DXY
all_dxy_fil$pop <- paste(all_dxy_fil$pop1, all_dxy_fil$pop2, sep = "-")
all_dxy_fil$pop_clean <- recode_factor(all_dxy_fil$pop, "glob_MR-glob_pure" = "MR globulus - Pure globulus", "glob_MR-cord_MR" = "MR globulus - cordata", "glob_pure-cord_MR" = "Pure globulus - cordata")
dxy_means <- data.frame(pop_clean = c("MR globulus - Pure globulus", "MR globulus - cordata", "Pure globulus - cordata"), mean = c(ref_mr_dxy, mr_cord_dxy, ref_cord_dxy))
dxy_plot <- ggplot(all_dxy_fil, aes(x = avg_dxy, fill = pop_clean)) + geom_histogram(color = "black") + 
           scale_fill_manual(name = "", values = c("#908168", "#45767e", "#3fe67c")) +
           facet_grid(pop_clean ~ .) + theme_bw() + labs(x = "Average Dxy per 5kb Window", y = "Count") +
           geom_vline(data=dxy_means, aes(xintercept = mean), colour="black") + 
           theme(legend.position = "none", axis.title=element_text(size=14))

# FST
all_fst_fil$pop <- paste(all_fst_fil$pop1, all_fst_fil$pop2, sep = "-")
all_fst_fil$pop_clean <- recode_factor(all_fst_fil$pop, "glob_MR-glob_pure" = "MR globulus - Pure globulus", "glob_MR-cord_MR" = "MR globulus - cordata", "glob_pure-cord_MR" = "Pure globulus - cordata")
fst_means <- data.frame(pop_clean = c("MR globulus - Pure globulus", "MR globulus - cordata", "Pure globulus - cordata"), mean = c(ref_mr_fst, mr_cord_fst, ref_cord_fst))
fst_plot <- ggplot(all_fst_fil, aes(x = avg_wc_fst, fill = pop_clean)) + geom_histogram(color = "black") + 
           scale_fill_manual(name = "", values = c("#908168", "#45767e", "#3fe67c")) +
           facet_grid(pop_clean ~ .) + theme_bw() + labs(x = "Average FST per 5kb Window", y = "Count") +
           geom_vline(data=fst_means, aes(xintercept = mean), colour="black") + 
           theme(legend.position = "none", axis.title=element_text(size=14))
```

## Tajima's D and fixed differences
Calculated Tajima's D in 50,000 bp sliding windows using `vcftools`.

```bash
# Run on UFRC queue system; see tajd_windows.job for more details.
# Resources used: 9 Mb, 4 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/"

vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/Ecordata.txt --TajimaD 50000 --out cord
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/Eglobulus_MR.txt --TajimaD 50000 --out globMR
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/Eglobulus_ref.txt --TajimaD 50000 --out globref
```

Calculated pi for individual sites for both reference populations using `vcftools`. I was looking for individual sites that have pi = 0 within species but pi != 0 between species as a signal of fixed differences. This excludes multi-allelic sites where one species has an allele fixed and the other has multiple other alleles unfixed, but given that biallelic SNPs make up ~1 million of the ~1.2 million SNPs total, I don't expect that particular scenario to seriously bias my results.

```bash
# Run on UFRC queue system; see pi_sites.job for more details.
# Resources used: 10 Mb, 11 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/Ecordata.txt --site-pi --out cord
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/Eglobulus_ref.txt --site-pi --out globref
vcftools --gzvcf "$VCF_DIR"/all_fil.vcf.gz --keep "$LIST_DIR"/ref_samples.txt --site-pi --out refs
```

Identified sites where pi = 0 within species and pi != 0 between species using `R`.
```bash
# Run on UFRC queue system; see fdiffs.job for more details.
# Resources used: 2.5 Gb, 2 days
module load R/4.3
Rscript get_fdiffs.r
```

See `chr_plots` README for synthesis of the results and plotting notes.