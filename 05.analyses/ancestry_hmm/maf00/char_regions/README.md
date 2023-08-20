# Characterize intervals of interest from Ancestry_HMM

## Delimit Introgressed Regions

Extracted introgressed regions identified by `Ancestry_HMM` for each sample using a custom `R` script. Introgressed regions were defined as:
* Containing no more than 750 bp in a row with less than 0.90 posterior probability of homozygous _E. cordata_ ancestry
* Containing at least one variant with a > 0.95 posterior probability of homozygous _E. cordata_ ancestry

```bash
# Performed in UFRC queue system; see get_intr_regions.job for more details.
# Resources used: 1.59 Mb, 1 min

module load R/4.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
POST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/posteriors"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    echo doing "$NAME"
    Rscript "$SCRIPT_DIR"/get_intr_bed.r "$POST_DIR"/"$NAME".posterior 0.95 0.90 750 0.85 ind_region_beds/"$NAME"_ahmm_regions.bed
done < "$LIST_DIR"/Eglobulus_MR.txt
```

For the regions on Chromosomes 6, 7, and 8 with many individuals sharing regions, retrieved shared regions in BED files by hand. Compared against heatmap plots and selected delimitation based on the blocks shared by most individuals; to resolve any small differences, took coordinate which would result in the smallest region.

| Chromosome | Start       | End           | Start taken from | End taken from |
| ---------- | ----------- | ------------- | ---------------- | -------------- |
| Chr06      | 13165829    | 13179404      | WC02/4229        | WC02/4229      |
| Chr07      | 234943      | 444947        | WA01/4190        | WA01/4190      |
| Chr07      | 7094890     | 7112122       | WG03/4505        | WG03/4505      |
| Chr08      | 22467118    | 22483317      | WG05/4241        | WG05/4241      |

Included the larger region for the first introgressed region in Chromosome 7 in the shared set because all but two samples had it.

## SNP Calling Stats
Checked variant info for marked regions as a quality check of introgression identification.

```bash
# Performed in UFRC's queue system; see call_stats.job for more details.
# Resources used: 5 Mb, 10 min
module load vcftools/0.1.16 
BEDDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/char_regions"
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

# Chromosome 6
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr06 --remove-indv SRR10339635 --bed "$BEDDIR"/chr06_shared.bed --site-mean-depth
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr06 --remove-indv SRR10339635 --bed "$BEDDIR"/chr06_shared.bed --site-quality
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr06 --remove-indv SRR10339635 --bed "$BEDDIR"/chr06_shared.bed --missing-site
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr06 --remove-indv SRR10339635 --bed "$BEDDIR"/chr06_shared.bed --get-INFO GL

# Chromosome 7
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr07 --remove-indv SRR10339635 --bed "$BEDDIR"/chr07_shared.bed --site-mean-depth
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr07 --remove-indv SRR10339635 --bed "$BEDDIR"/chr07_shared.bed --site-quality
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr07 --remove-indv SRR10339635 --bed "$BEDDIR"/chr07_shared.bed --missing-site
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr07 --remove-indv SRR10339635 --bed "$BEDDIR"/chr07_shared.bed --get-INFO GL

# Chromosome 8
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr08 --remove-indv SRR10339635 --bed "$BEDDIR"/chr08_shared.bed --site-mean-depth
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr08 --remove-indv SRR10339635 --bed "$BEDDIR"/chr08_shared.bed --site-quality
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr08 --remove-indv SRR10339635 --bed "$BEDDIR"/chr08_shared.bed --missing-site
vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --out all_ahmm_chr08 --remove-indv SRR10339635 --bed "$BEDDIR"/chr08_shared.bed --get-INFO GL
```

Plotted results in `R`.

```R
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/char_regions/call_stats"

setwd(wdir)
# Depth
chr06_dep_tab <- read.table("all_ahmm_chr06.ldepth.mean", sep = "\t", header = TRUE)
plot(chr06_dep_tab$POS, chr06_dep_tab$MEAN_DEPTH, type = "l", main = "Chr06 Average Depth", ylab = "Avg Depth per sample (reads)", xlab = "Position")
abline(h = mean(chr06_dep_tab$MEAN_DEPTH), col = "red")
text(x = 13178500, y = 38.75, paste("Mean:", round(mean(chr06_dep_tab$MEAN_DEPTH), 3)), col = "red")

chr07_dep_tab <- read.table("all_ahmm_chr07.ldepth.mean", sep = "\t", header = TRUE)
chr07_1 <- which(chr07_dep_tab$POS < 500000)
plot(chr07_dep_tab$POS[chr07_1], chr07_dep_tab$MEAN_DEPTH[chr07_1], type = "l", main = "Chr07-1 Average Depth", ylab = "Avg Depth per sample (reads)", xlab = "Position")
abline(h = mean(chr07_dep_tab$MEAN_DEPTH[chr07_1]), col = "red")
text(x = 400000, y = 47, paste("Mean:", round(mean(chr07_dep_tab$MEAN_DEPTH[chr07_1]), 3)), col = "red")

chr07_2 <- which(chr07_dep_tab$POS > 500000)
plot(chr07_dep_tab$POS[chr07_2], chr07_dep_tab$MEAN_DEPTH[chr07_2], type = "l", main = "Chr07-2 Average Depth", ylab = "Avg Depth per sample (reads)", xlab = "Position")
abline(h = mean(chr07_dep_tab$MEAN_DEPTH[chr07_2]), col = "red")
text(x = 7110000, y = 35, paste("Mean:", round(mean(chr07_dep_tab$MEAN_DEPTH[chr07_2]), 3)), col = "red")

chr08_dep_tab <- read.table("all_ahmm_chr08.ldepth.mean", sep = "\t", header = TRUE)
plot(chr08_dep_tab$POS, chr08_dep_tab$MEAN_DEPTH, type = "l", main = "Chr08 Average Depth", ylab = "Avg Depth per sample (reads)", xlab = "Position")
abline(h = mean(chr08_dep_tab$MEAN_DEPTH), col = "red")
text(x = 22470000, y = 43, paste("Mean:", round(mean(chr08_dep_tab$MEAN_DEPTH), 3)), col = "red")

# % Missing
chr06_mis_tab <- read.table("all_ahmm_chr06.lmiss", sep = "\t", header = TRUE)
plot(chr06_mis_tab$POS, chr06_mis_tab$F_MISS, type = "l", main = "Chr06 Fraction of GTs Missing", ylab = "Fraction Missing", xlab = "Position")
abline(h = mean(chr06_mis_tab$F_MISS), col = "red")
text(x = 13178000, y = 0.12, paste("Mean:", round(mean(chr06_mis_tab$F_MISS, na.rm = TRUE), 3)), col = "red")

chr07_mis_tab <- read.table("all_ahmm_chr07.lmiss", sep = "\t", header = TRUE)
chr07_1 <- which(chr07_mis_tab$POS < 500000)
plot(chr07_mis_tab$POS[chr07_1], chr07_mis_tab$F_MISS[chr07_1], type = "l", main = "Chr07-1 Fraction of GTs Missing", ylab = "Fraction Missing", xlab = "Position")
abline(h = mean(chr07_mis_tab$F_MISS[chr07_1]), col = "red")
text(x = 400000, y = 0.12, paste("Mean:", round(mean(chr07_mis_tab$F_MISS[chr07_1]), 3)), col = "red")

chr07_2 <- which(chr07_mis_tab$POS > 500000)
plot(chr07_mis_tab$POS[chr07_2], chr07_mis_tab$F_MISS[chr07_2], type = "l", main = "Chr07-2 Fraction of GTs Missing", ylab = "Fraction Missing", xlab = "Position")
abline(h = mean(chr07_mis_tab$F_MISS[chr07_2]), col = "red")
text(x = 7110000, y = 0.12, paste("Mean:", round(mean(chr07_mis_tab$F_MISS[chr07_2]), 3)), col = "red")

chr08_mis_tab <- read.table("all_ahmm_chr08.lmiss", sep = "\t", header = TRUE)
plot(chr08_mis_tab$POS, chr08_mis_tab$F_MISS, type = "l", main = "Chr08 Fraction GTs Missing", ylab = "Fraction Missing", xlab = "Position")
abline(h = mean(chr08_mis_tab$F_MISS), col = "red")
text(x = 22470000, y = 0.12, paste("Mean:", round(mean(chr08_mis_tab$F_MISS), 3)), col = "red")

# QUAL
chr06_qual_tab <- read.table("all_ahmm_chr06.lqual", sep = "\t", header = TRUE)
plot(chr06_qual_tab$POS, chr06_qual_tab$QUAL, type = "l", main = "Chr06 Quality", ylab = "Quality Score", xlab = "Position")
abline(h = mean(chr06_qual_tab$QUAL), col = "red")
text(x = 13176000, y = 500, paste("Mean:", round(mean(chr06_qual_tab$QUAL, na.rm = TRUE), 0)), col = "red")

chr07_qual_tab <- read.table("all_ahmm_chr07.lqual", sep = "\t", header = TRUE)
chr07_1 <- which(chr07_qual_tab$POS < 500000)
plot(chr07_qual_tab$POS[chr07_1], chr07_qual_tab$QUAL[chr07_1], type = "l", main = "Chr07-1 Quality", ylab = "Quality Score", xlab = "Position")
abline(h = mean(chr07_qual_tab$QUAL[chr07_1]), col = "red")
text(x = 400000, y = 500, paste("Mean:", round(mean(chr07_qual_tab$QUAL[chr07_1]), 0)), col = "red")

chr07_2 <- which(chr07_mis_tab$POS > 500000)
plot(chr07_qual_tab$POS[chr07_2], chr07_qual_tab$QUAL[chr07_2], type = "l", main = "Chr07-2 Quality", ylab = "Quality Score", xlab = "Position")
abline(h = mean(chr07_qual_tab$QUAL[chr07_2]), col = "red")
text(x = 7110000, y = 500, paste("Mean:", round(mean(chr07_qual_tab$QUAL[chr07_2]), 0)), col = "red")

chr08_qual_tab <- read.table("all_ahmm_chr08.lqual", sep = "\t", header = TRUE)
plot(chr08_qual_tab$POS, chr08_qual_tab$QUAL, type = "l", main = "Chr08 Quality", ylab = "Quality Score", xlab = "Position")
abline(h = mean(chr08_qual_tab$QUAL), col = "red")
text(x = 22470000, y = 500, paste("Mean:", round(mean(chr08_qual_tab$QUAL), 0)), col = "red")
```

## genome scan statistics
Checked pi, dxy, and f-stats for intervals.

```R
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/char_regions/gscan_stats"
dxy_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy/all_dxy.txt"
fst_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy/all_fst.txt"
pi_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy/all_pi.txt"
fstats_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/dsuite/localFstats_40_20.txt"

setwd(wdir)
dxy_tab <- read.table(dxy_filename, header = TRUE, na.strings = c("NA"), sep = "\t")
fst_tab <- read.table(fst_filename, header = TRUE, na.strings = c("NA"), sep = "\t")
pi_tab <- read.table(pi_filename, header = TRUE, na.strings = c("NA"), sep = "\t")
fstats_tab <- read.table(fstats_filename, header = TRUE, sep = "\t")

ahmm_regions <- data.frame(chr = c("Chr06", "Chr07", "Chr07", "Chr08"), start = c(13165829, 234943, 7094890, 22467118), end = c(13179404, 444947, 7112122, 22483317))

# pi
pi_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(pi_tab$chromosome == ahmm_regions[i, "chr"])
    pop_rows <- which(pi_tab$pop == "glob_MR")
    good_cov <- which(pi_tab$no_sites > 40)
    starts <- which((pi_tab$window_pos_1 > ahmm_regions[i, "start"]) & (pi_tab$window_pos_1 < ahmm_regions[i, "end"]))
    ends <- which((pi_tab$window_pos_2 > ahmm_regions[i, "start"]) & (pi_tab$window_pos_2 < ahmm_regions[i, "end"]))
    in_region <- sort(unique(union(starts, ends)))
    rows_to_calc <- intersect(in_region, intersect(good_cov, intersect(chr_rows, pop_rows)))

    pi_val <- sum(pi_tab[rows_to_calc, "count_diffs"])/sum(pi_tab[rows_to_calc, "count_comparisons"])
    pi_vec <- c(pi_vec, pi_val)
}

# dxy (cordata)
dxy_cord_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(dxy_tab$chromosome == ahmm_regions[i, "chr"])
    pop_rows <- which((dxy_tab$pop1 == "glob_MR") & (dxy_tab$pop2 == "cord_MR"))
    good_cov <- which(dxy_tab$no_sites > 40)
    starts <- which((dxy_tab$window_pos_1 > ahmm_regions[i, "start"]) & (dxy_tab$window_pos_1 < ahmm_regions[i, "end"]))
    ends <- which((dxy_tab$window_pos_2 > ahmm_regions[i, "start"]) & (dxy_tab$window_pos_2 < ahmm_regions[i, "end"]))
    in_region <- sort(unique(union(starts, ends)))
    rows_to_calc <- intersect(in_region, intersect(good_cov, intersect(chr_rows, pop_rows)))

    dxy_cord_val <- sum(dxy_tab[rows_to_calc, "count_diffs"])/sum(dxy_tab[rows_to_calc, "count_comparisons"])
    dxy_cord_vec <- c(dxy_cord_vec, dxy_cord_val)
}

# dxy (pure globulus)
dxy_glob_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(dxy_tab$chromosome == ahmm_regions[i, "chr"])
    pop_rows <- which((dxy_tab$pop1 == "glob_MR") & (dxy_tab$pop2 == "glob_pure"))
    good_cov <- which(dxy_tab$no_sites > 40)
    starts <- which((dxy_tab$window_pos_1 > ahmm_regions[i, "start"]) & (dxy_tab$window_pos_1 < ahmm_regions[i, "end"]))
    ends <- which((dxy_tab$window_pos_2 > ahmm_regions[i, "start"]) & (dxy_tab$window_pos_2 < ahmm_regions[i, "end"]))
    in_region <- sort(unique(union(starts, ends)))
    rows_to_calc <- intersect(in_region, intersect(good_cov, intersect(chr_rows, pop_rows)))

    dxy_glob_val <- sum(dxy_tab[rows_to_calc, "count_diffs"])/sum(dxy_tab[rows_to_calc, "count_comparisons"])
    dxy_glob_vec <- c(dxy_glob_vec, dxy_glob_val)
}

# fst (cord)
fst_cord_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(fst_tab$chromosome == ahmm_regions[i, "chr"])
    pop_rows <- which((fst_tab$pop1 == "glob_MR") & (fst_tab$pop2 == "cord_MR"))
    good_cov <- which(fst_tab$no_snps > 40)
    starts <- which((fst_tab$window_pos_1 > ahmm_regions[i, "start"]) & (fst_tab$window_pos_1 < ahmm_regions[i, "end"]))
    ends <- which((fst_tab$window_pos_2 > ahmm_regions[i, "start"]) & (fst_tab$window_pos_2 < ahmm_regions[i, "end"]))
    in_region <- sort(unique(union(starts, ends)))
    rows_to_calc <- intersect(in_region, intersect(good_cov, intersect(chr_rows, pop_rows)))

    fst_cord_val <- mean(fst_tab[rows_to_calc, "avg_wc_fst"], na.rm = TRUE)
    fst_cord_vec <- c(fst_cord_vec, fst_cord_val)
}

# fst (glob)
fst_glob_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(fst_tab$chromosome == ahmm_regions[i, "chr"])
    pop_rows <- which((fst_tab$pop1 == "glob_MR") & (fst_tab$pop2 == "glob_pure"))
    good_cov <- which(fst_tab$no_snps > 40)
    starts <- which((fst_tab$window_pos_1 > ahmm_regions[i, "start"]) & (fst_tab$window_pos_1 < ahmm_regions[i, "end"]))
    ends <- which((fst_tab$window_pos_2 > ahmm_regions[i, "start"]) & (fst_tab$window_pos_2 < ahmm_regions[i, "end"]))
    in_region <- sort(unique(union(starts, ends)))
    rows_to_calc <- intersect(in_region, intersect(good_cov, intersect(chr_rows, pop_rows)))

    fst_glob_val <- mean(fst_tab[rows_to_calc, "avg_wc_fst"], na.rm = TRUE)
    fst_glob_vec <- c(fst_glob_vec, fst_glob_val)
}


# f-stats
fdm_vec <- c()
df_vec <- c()
for(i in c(1:nrow(ahmm_regions))){
    chr_rows <- which(fstats_tab$chr == ahmm_regions[i, "chr"])
    starts <- which(fstats_tab$windowStart < ahmm_regions[i, "start"])
    ends <- which(fstats_tab$windowEnd > ahmm_regions[i, "end"])
    in_region <- sort(unique(intersect(starts, ends)))
    rows_to_calc <- intersect(chr_rows, in_region)

    fdm_val <- mean(fstats_tab[rows_to_calc, "f_dM"], na.rm = TRUE)
    fdm_vec <- c(fdm_vec, fdm_val)
    df_val <- mean(fstats_tab[rows_to_calc, "d_f"], na.rm = TRUE)
    df_vec <- c(df_vec, df_val)
}
```

| Chromosome | Start (bp)  | End (bp)    | pi          | dxy (cord)   | dxy (glob pure) | fDm          | df           | fst (cord) | fst (glob) |
| ---------- | ----------- | ----------- | ----------- | ------------ | --------------- | ------------ | ------------ | ---------- | ---------- |
| Chr06      | 13165829    | 13179404    | 0.05609     | 0.04972      | 0.07441         | 0.2167       | 0.09257      | 0.1126     | 0.06604    |
| Chr07      | 234943      | 444947      | 0.06381     | 0.05954      | 0.06996         | NA           | NA           | 0.1542     | 0.03204    |
| Chr07      | 7094890     | 7112122     | 0.03139     | 0.02401      | 0.03287         | 0.08769      | 0.01938      | 0.1268     | -0.01993   |
| Chr08      | 22467118    | 22483317    | 0.1272      | 0.1576       | 0.1128          | 0.01257      | -0.04920     | NA         | NA         |

All segments have some statistic where values make sense under introgression, but other statistics which do not, though I haven't tested against the actual distributions of each statistic. There are several regions near the edge of the Chromosome 7 arm that did not have enough information to calculate pi/dxy, which is somewhat concerning and needs to be looked into.

Plotted distribution of pi, dxy, fDm, df and average value of those statistics within introgressed windows.

```R

```