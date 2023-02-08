# Whole-genome PCA
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/pca/).

## Prune linked SNPs
Note: this step will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see link_prune.job for details
# Resources used: 2 Mb, 30 sec

module load plink/1.90b3.39 

export INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
export INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"

plink --vcf "$INFILE00" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --vcf-half-call m --out maf0.00/all_maf0.00
plink --vcf "$INFILE05" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --vcf-half-call m --out maf0.05/all_maf0.05

```

## Calculate principal components
Note: this step also will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see pca.job for details
# Resources used: 2 Mb, 30 sec

module load plink/1.90b3.39 

export INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
export INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"

plink --vcf "$INFILE00" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract all_maf0.00.prune.in --vcf-half-call m --make-bed --pca --out all_maf0.00
plink --vcf "$INFILE05" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract all_maf0.05.prune.in --vcf-half-call m --make-bed --pca --out all_maf0.05
```

## Plot principal components
Performed locally in [`R programming language`](https://www.r-project.org/). Used metadata to color samples by type: reference _E. cordata_ ("cord_MR", yellow), reference _E. globulus_ ("glob_ref", blue), or introgressed _E. globulus_ ("glob_MR", black).

**MAF=0.00:**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00/all_maf0.00.eigenvec"
eigenval_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00/all_maf0.00.eigenval"

setwd(working_dir)

# Read input files
sample_table <- read.csv(table_address, header = TRUE, as.is = TRUE)
pca <- read.table(pca_address, header = FALSE)
eigenval <- read.table(eigenval_address, header = FALSE)$V1

# Reformat imported data
pca <- pca[,-2]
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste("PC", c(1:(ncol(pca)-1)), sep = "")
pca <- cbind(pca, sample_table[match(pca$sample, sample_table$RAPiD_ID), "Taxon"])
names(pca)[ncol(pca)] <- "taxon"

# Plot eigenvalues
pve <- data.frame(PC = 1:(ncol(pca) - 2), pve = eigenval/sum(eigenval))
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00, PC1 vs. PC2")
b
# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8) + geom_text(data = subset(pca, sample == "WG04"), aes(label = sample), hjust = -0.25, vjust = -0.25)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00, PC1 vs. PC3")
c
# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + geom_text(data = subset(pca, sample %in% c("WF03", "WG04")), aes(label = sample), hjust = 1, vjust = -.25)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00, PC2 vs. PC3")
d
```

**MAF=0.05:**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.05"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.05/all_maf0.05.eigenvec"
eigenval_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.05/all_maf0.05.eigenval"

setwd(working_dir)

# Read input files
sample_table <- read.csv(table_address, header = TRUE, as.is = TRUE)
pca <- read.table(pca_address, header = FALSE)
eigenval <- read.table(eigenval_address, header = FALSE)$V1

# Reformat imported data
pca <- pca[,-2]
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste("PC", c(1:(ncol(pca)-1)), sep = "")
pca <- cbind(pca, sample_table[match(pca$sample, sample_table$RAPiD_ID), "Taxon"])
names(pca)[ncol(pca)] <- "taxon"

# Plot eigenvalues
pve <- data.frame(PC = 1:(ncol(pca) - 2), pve = eigenval/sum(eigenval))
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.05")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.6) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05, PC1 vs. PC2")
b
# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05, PC1 vs. PC3")
c
# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + xlim(-0.75, 0.35) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05, PC2 vs. PC3")
d
```

**Results:**

Two samples, WF03/1051 and WG04/2025, were very different from the rest of their groups (reference _E. globulus_ and introgressed _E. globulus_ respectively) for MAF=0.00 plots. WG04/2025 fell back within range of the rest of introgressed _E. globulus_ variation in MAF=0.05 plots. Below, MAF=0.05 plots are shown.


![Plot of PC1 vs PC2; PC1 differentiates strongly between _E. cordata_ and _E. globulus_ (both introgressant and pure), while PC2 differentiates between introgressant and pure _E. globulus_, with _E. cordata_ clustering with introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc12.png "PC1 vs. PC2")

![Plot of PC1 vs PC3; very similar to PC1 vs PC2 but _E. globulus_ populations intergrade more along PC3.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc13.png "PC1 vs. PC3")

![Plot of PC2 vs PC3; pure _E. globulus_ seems to cluster slightly away from the other individuals, while _E. cordata_ clusters tightly in the center of a cloud of introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc23.png "PC2 vs. PC3")

![Barplot of percentage genetic variation explained by each principal component; PC1 explains about 30% of variation, while all other PCs pictured explain about 5%](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc_var_explained.png "Percent Variance Explained by each PC")

## Investigate Odd Samples
WF03/1051 and WG04/2025 were distant in 2nd and 3rd PC plots from the rest of their sample groups. Calculated Fsts from odd samples to the rest of the groups to decide how big of an issue this would be.

```bash
module load vcftools/0.1.16

INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
LISTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/weird_sample_check"

# WF03, MAF=0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_wf03_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.00_wf03_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.00_wf03_cord

# WG04, MAF=0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_wg04_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --out maf0.00_wg04_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.00_wg04_cord

# WF03, MAF=0.05
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_wf03_globref
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.05_wf03_globmr
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.05_wf03_cord
```

Get genome-wide average in `R`:

```R
maf00_wf03_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/weird_sample_check/maf0.00_wf03"
maf00_wg04_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/weird_sample_check/maf0.00_wg04"
maf05_wf03_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/weird_sample_check/maf0.05_wf03"

# MAF=0.00, WF03 versus the three sample groups
setwd(maf00_wf03_dir)
maf00_wf03_globref_tab <- read.table("maf0.00_wf03_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_wf03_globmr_tab <- read.table("maf0.00_wf03_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_wf03_cord_tab <- read.table("maf0.00_wf03_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_cord_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.00, WG04 versus the three sample groups
setwd(maf00_wg04_dir)
maf00_wg04_globref_tab <- read.table("maf0.00_wg04_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_wg04_globmr_tab <- read.table("maf0.00_wg04_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_wg04_cord_tab <- read.table("maf0.00_wg04_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_cord_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, WF03 versus the three sample groups
setwd(maf05_wf03_dir)
maf05_wf03_globref_tab <- read.table("maf0.05_wf03_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_globref_tab$WEIR_AND_COCKERHAM_FST))
maf05_wf03_globmr_tab <- read.table("maf0.05_wf03_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf05_wf03_cord_tab <- read.table("maf0.05_wf03_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_cord_tab$WEIR_AND_COCKERHAM_FST))
```

Results:

| MAF  | Group 1   | Group 2      | Min   | 1st Quart | Median | Mean  | 3rd Quart | Max  |
|------|-----------|--------------|-------|-----------|--------|------ |-----------|------|
| 0.00 | WF03/1051 | E. glob ref  | -1.09 | -0.36     | -0.22  | -0.11 | 0.056     | 1.00 |
| 0.00 | WF03/1051 | E. glob intr | -1.00 | -0.31     | -0.22  | -0.11 | 0.018     | 1.00 |
| 0.00 | WF03/1051 | E. cordata   | -1.09 | -0.20     | 0.051  | 0.19  | 0.72      | 1.00 |
| 0.00 | WG04/2025 | E. glob ref  | -1.09 | -0.36     | -0.26  | -0.13 | 0.00      | 1.00 |
| 0.00 | WG04/2025 | E. glob intr | -1.00 | -0.34     | -0.23  | -0.13 | -0.03     | 1.00 |
| 0.00 | WG04/2025 | E. cordata   | -1.19 | -0.22     | 0.052  | 0.19  | 0.76      | 1.00 |
| 0.05 | WF03/1051 | E. glob ref  | -1.09 | -0.36     | -0.18  | -0.094| 0.074     | 1.00 |
| 0.05 | WF03/1051 | E. glob intr | -1.00 | -0.28     | -0.16  | -0.077| 0.062     | 1.00 |
| 0.05 | WF03/1051 | E. cordata   | -1.09 | -0.099    | 0.17   | 0.27  | 0.82      | 1.00 |

**Conclusion:** Divergence is pretty consistent between problem samples, their own group, and the other _E. globulus_ group but is much less than for problem samples with _E. cordata_. I can't see a reason to exclude them since they're exhibiting the pattern of divergence I'd expect of any other sample (can't speak to the degree of divergence though). HOWEVER: Since most of Fst values are negative (indicating larger variation within than between populations), it might not be the best metric to use for evaluating these samples. Will need to test whether results are robust to excluding WF03/1051 at the very least.
