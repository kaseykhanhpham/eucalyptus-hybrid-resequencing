# PC-correlated SNPs

Identified SNPs in reference populations excessively associated with genome-wide principle component for **MAF=0.05** using the `R` package [`pcadapt`](https://cran.r-project.org/web/packages/pcadapt/index.html).

Referred to [this tutorial](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html) from the developers.

## Get correlated SNPs

Subsetted all-chromosome just-variants VCF to reference samples only using `bcftools`. ~~Pruned linked variants~~, converted to BED file, and calculated PCs using `PLINK`.

```bash
module load bcftools/1.15
module load plink/1.90b3.39
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_var"

ls "$VCFDIR"/*.vcf.gz > files_to_concat.txt
# bcftools concat --file-list files_to_concat.txt --output-type z | bcftools view --samples-file ref_samples.txt --min-af 0.05 --output-type z --output "$VCFDIR"/refs_fil_maf05.vcf.gz

bcftools concat --file-list files_to_concat.txt --output-type z | bcftools view --samples-file ref_samples.txt --min-af 0.00 --output-type z --output "$VCFDIR"/refs_fil_maf00.vcf.gz

# plink --vcf "$VCFDIR"/refs_fil_maf05.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --vcf-half-call m --out refs_maf05
# plink --vcf "$VCFDIR"/refs_fil_maf05.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract refs_maf05.prune.in --vcf-half-call m --make-bed --pca --out refs_maf05

plink --vcf "$VCFDIR"/refs_fil_maf00.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call m --make-bed --pca --out refs_maf00
```

### Decide on PCs to analyze

Decided on K value to use based on plots of PCs. PCs which differentiated between reference groups were selected.

```R
# Performed on local computer.
library(ggplot2)
# set workspace vars
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/pca_corr"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
# pca_address <- paste(wdir, "refs_maf05.eigenvec", sep = "/")
# eigenval_address <- paste(wdir, "refs_maf05.eigenval", sep = "/")
pca_address <- paste(wdir, "refs_maf00.eigenvec", sep = "/")
eigenval_address <- paste(wdir, "refs_maf00.eigenval", sep = "/")
setwd(wdir)

# import files
sample_table <- read.csv(table_address, header = TRUE, as.is = TRUE)
pca <- read.table(pca_address, header = FALSE)
eigenval <- read.table(eigenval_address, header = FALSE)$V1

# format PCA table and add a column for taxon identifiers
pca <- pca[,-2]
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste("PC", c(1:(ncol(pca)-1)), sep = "")
pca <- cbind(pca, sample_table[match(pca$sample, sample_table$RAPiD_ID), "Taxon"])
names(pca)[ncol(pca)] <- "taxon"

# Plot eigenvalues
pve <- data.frame(PC = 1:(ncol(pca) - 2), pve = eigenval/sum(eigenval))
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
# a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.05, reference populations only")
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00, reference populations only")
cumsum(pve$pve)

# Plot PCs
## PC1 vs. PC2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + scale_colour_manual(values = c("goldenrod1", "deepskyblue4")) + theme_light() + coord_equal()
# b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 just ref pops, PC1 vs. PC2")
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 just ref pops, PC1 vs. PC2")
b
## PC3 vs. PC4
c <- ggplot(pca, aes(PC3, PC4, col = taxon)) + geom_point(size = 2) + scale_colour_manual(values = c("goldenrod1", "deepskyblue4")) + theme_light() + coord_equal()
# c <- c + xlab(paste("PC3 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC4 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 just ref pops, PC3 vs. PC4")
c <- c + xlab(paste("PC3 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC4 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 just ref pops, PC3 vs. PC4")
c
# PC5 vs. PC6
d <- ggplot(pca, aes(PC5, PC6, col = taxon)) + geom_point(size = 2) + scale_colour_manual(values = c("goldenrod1", "deepskyblue4")) + theme_light() + coord_equal()
# d <- d + xlab(paste("PC5 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC6 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 just ref pops, PC5 vs. PC6")
d <- d + xlab(paste("PC5 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC6 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 just ref pops, PC5 vs. PC6")
d
```

PC1 is the only PC which strongly differentiates between _E. globulus_ and _E. cordata_, just as in all-sample analysis. Therefore, will only look for variants correlated with PC1 (K = 1).

![Plot of PC1 vs PC2 for MAF=0.00 (just reference populations); PC1 differentiates strongly between _E. cordata_ and _E. globulus_, while _E. cordata_ clusters in the middle of a wider _E. globulus_ distribution along PC2.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/pca_corr/refs_maf00_pc12.png "MAF=0.00, refs only, PC1 vs. PC2")

![Plot of PC3 vs PC4 for MAF=0.00 (just reference populations); _E. cordata_ clusters in the middle of a wider _E. globulus_ distribution along PCs 3 and 4.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/pca_corr/refs_maf00_pc34.png "MAF=0.00, refs only, PC3 vs. PC4")

![Plot of PC5 vs PC6 for MAF=0.00 (just reference populations); _E. cordata_ clusters in the middle of a wider _E. globulus_ distribution along PCs 5 and 6.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/pca_corr/refs_maf00_pc56.png "MAF=0.00, refs only, PC5 vs. PC6")

### Get outliers

```R
library(pcadapt)
library(qvalue)
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/pca_corr"
# bim_address <- paste(wdir, "refs_maf05.bim", sep = "/")
bim_address <- paste(wdir, "refs_maf00.bim", sep = "/")

setwd(wdir)
# bed_in <- read.pcadapt(input = "refs_maf05.bed", type = "bed")
bed_in <- read.pcadapt(input = "refs_maf00.bed", type = "bed")
bim_in <- read.table(bim_address, sep = "\t", header = FALSE)
corr_test <- pcadapt(input = bed_in, K = 1, min.maf = 0.000001)
summary(corr_test$pvalues)
plot(corr_test, option = "manhattan")
plot(corr_test, option = "qqplot")
hist(corr_test$pvalues, xlab = "p-values", main = NULL, breaks = 50)
plot(corr_test, option = "stat.distribution")

## Get outliers
qval <- qvalue(corr_test$pvalues)$qvalues
alpha <- 0.001
outliers <- which(qval < alpha)

bim_out <- bim_in[outliers,]
# write.table(bim_out, "refs_maf05_outliers.bim", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(bim_out, "refs_maf00_outliers.bim", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
```

There were 22,532 SNPs out of 1,145,219 correlated with PC1 after adjustment.

Checked PC-correlated SNPs against outlier windows identified in genome scan analysis.

```bash
module load R/4.2
module load bedtools/2.30.0
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
GS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

# Rscript "$SCRIPT_DIR"/bim_to_bed.r refs_maf05_outliers.bim refs_maf05_outliers.bed
Rscript "$SCRIPT_DIR"/bim_to_bed.r refs_maf00_outliers.bim refs_maf00_outliers.bed
# bedtools intersect -a "$GS_DIR"/pi_dxy_outl_p15.bed -b refs_maf05_outliers.bed -wb > pi15_dxy15_pca_overlap.bed
bedtools intersect -a "$GS_DIR"/pi_dxy_outl_p15.bed -b refs_maf00_outliers.bed -wb > pi15_dxy15_pca_overlap.bed
```

No PC-correlated variants overlapping with pi/dxy outlier windows.

## Compare PC-correlated SNPs in introgressed samples
