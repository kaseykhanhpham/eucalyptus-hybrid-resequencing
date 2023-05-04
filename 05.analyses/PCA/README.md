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

**MAF=0.00, without outgroup**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00/all_maf0.00_ingr.eigenvec"
eigenval_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00/all_maf0.00_ingr.eigenval"

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
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 just ingr, PC1 vs. PC2")
b
b + xlim(-0.125, 0.30) + ylim(-0.125, 0.3) + ggtitle("maf=0.00 just ingr, PC1 vs. PC2, zoomed") # zoomed to ignore outliers

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8) + geom_text(data = subset(pca, sample == "WG04"), aes(label = sample), hjust = -0.25, vjust = -0.25)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 just ingr, PC1 vs. PC3")
c
c + xlim(-0.125, 0.3) + ylim(-0.125, 0.3) + ggtitle("maf=0.00 just ingr, PC1 vs. PC3, zoomed") # zoomed to ignore outliers

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + geom_text(data = subset(pca, sample %in% c("WF03", "WG04")), aes(label = sample), hjust = 1, vjust = -.25)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 just ingr, PC2 vs. PC3")
d
d + xlim(-0.075, 0.1) + ylim(-0.125, 0.05) + ggtitle("maf=0.00 just ingr, PC2 vs. PC3, zoomed") # zoomed to ignore outliers
```

**MAF=0.00, with outgroup:**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.00/with_outgroup"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf0.00.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf0.00.eigenval", sep = "/")

setwd(working_dir)

# Read input files
sample_table <- read.csv(table_address, header = TRUE, as.is = TRUE)
sample_table <- rbind(sample_table, c("SRR10339635", "outgroup", "NA"))
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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00 with outgroup")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + ylim(-0.4, 0.7)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 with outgr, PC1 vs. PC2")
b
b + xlim(-0.35, 0.10) + ylim(-0.125, 0.3) + ggtitle("maf=0.00 with outgr, PC1 vs. PC2, zoomed") # zoomed to ignore outgroup
b + xlim(-0.04, -0.02) + ylim(-0.105, -0.08) + geom_text(data = subset(pca, sample %in% c("WF03", "WG04")), aes(label = sample), hjust = 1, vjust = -.25) + ggtitle("maf=0.00 with outgr, PC1 vs. PC2, just E. globulus")

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 with outgr, PC1 vs. PC3")
c
c + xlim(-0.05, 0.05) + ylim(-0.05, 0.05) + ggtitle("maf=0.00 with outgr, PC1 vs. PC3, zoomed") # zoomed to ignore outgroup and outliers

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + xlim(-0.4, 0.7) + geom_point(size = 2) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = 1, vjust = -.25)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 with outgr, PC2 vs. PC3")
d
d + xlim(-0.125, 0.305) + ylim(-0.30, 0.1) + ggtitle("maf=0.00 with outgr, PC2 vs. PC3, zoomed") # zoomed to ignore outliers
```

**MAF=0.05, without outgroup:**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.05"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf0.05_ingr.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf0.05_ingr.eigenval", sep = "/")

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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.05 just ingr")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.7) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 just ingr, PC1 vs. PC2")
b
b + xlim(-0.2, 0.35) + ylim(-0.35, 0.20) + ggtitle("maf=0.05 just ingr, PC1 vs. PC2, zoomed") # zoomed to ignore outlier

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 just ingr, PC1 vs. PC3")
c
c + xlim(-0.125, 0.5) + ylim(-0.35, 0.3) + ggtitle("maf=0.05 just ingr, PC1 vs. PC3, zoomed") # zoomed to ignore outlier

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = 1, vjust = 1)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 just ingr, PC2 vs. PC3")
d
d + xlim(-0.375, 0.25) + ylim(-0.35, 0.3) + ggtitle("maf=0.05 just ingr, PC2 vs. PC3, zoomed")
```

**MAF=0.00, with outgroup:**

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/maf0.05/with_outgroup"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf0.05.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf0.05.eigenval", sep = "/")

setwd(working_dir)

# Read input files
sample_table <- read.csv(table_address, header = TRUE, as.is = TRUE)
sample_table <- rbind(sample_table, c("SRR10339635", "outgroup", "NA"))
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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00 with outgroup")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.4, 0.7)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 with outgr, PC1 vs. PC2")
b
b + xlim(-0.11, 0.3) + ylim(-0.11, 0.3) + ggtitle("maf=0.05 with outgr, PC1 vs. PC2, zoomed") # zoomed to ignore outgroup
b + xlim(-0.11, 0) + ylim(-0.05, 0.075) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = 1, vjust = -.25) + ggtitle("maf=0.00 with outgr, PC1 vs. PC2, just E. globulus")

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = -0.25, vjust = -0.25) + xlim(-0.45, 0.65)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 with outgr, PC1 vs. PC3")
c
c + xlim(-0.2, 0.3) + ylim(-0.2, 0.3) + ggtitle("maf=0.05 with outgr, PC1 vs. PC3, zoomed") # zoomed to ignore outliers

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + geom_text(data = subset(pca, sample == "WF03"), aes(label = sample), hjust = 1, vjust = -.25) + ylim(-0.7, 0.4)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4", "gray"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 with outgr, PC2 vs. PC3")
d
d + xlim(-0.125, 0.45) + ylim(-0.2, 0.375) + ggtitle("maf=0.05 with outgr, PC2 vs. PC3, zoomed") # zoomed to ignore outliers and outgroup
```

**Results:**

Two samples, WF03/1051 and WG04/2025, were very different from the rest of their groups (reference _E. globulus_ and introgressed _E. globulus_ respectively) for MAF=0.00 plots. WG04/2025 fell back within range of the rest of introgressed _E. globulus_ variation in MAF=0.05 plots. Otherwise, MAF=0.00 and MAF=0.05 plots were similar in pattern; samples were more tightly clustered for MAF=0.00, indicating that singletons generally tell the same story as shared variants. Inclusion of the outgroup also did not change the general pattern expressed, it only shunted the observed patterns to PCs 2 and 3, with distance between the ingroup and outgroup dominating PC1. Therefore, only MAF=0.05 plots excluding the outgroup are shown below. Plots with the outgroup and for MAF=0.00 can be viewed in their designated directories.
 

![Plot of PC1 vs PC2; PC1 differentiates strongly between _E. cordata_ and _E. globulus_ (both introgressant and pure), while PC2 differentiates between introgressant and pure _E. globulus_, with _E. cordata_ clustering with introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf0.05_ingr_pc12.png "PC1 vs. PC2")

![Plot of PC1 vs PC3; very similar to PC1 vs PC2 but _E. globulus_ populations intergrade more along PC3.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf0.05_ingr_pc13.png "PC1 vs. PC3")

![Plot of PC2 vs PC3; pure _E. globulus_ seems to cluster slightly away from the other individuals, while _E. cordata_ clusters tightly in the center of a cloud of introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf0.05_ingr_pc23.png "PC2 vs. PC3")

![Barplot of percentage genetic variation explained by each principal component; PC1 explains about 30% of variation, while all other PCs pictured explain about 5%](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf0.05_ingr_variance_explained.png "Percent Variance Explained by each PC")

## Outlier Sample Check
Checked PCA without outgroup and sample WF03.

### Calculate PCA eigenvalues for both MAF levels

```bash
# Performed on UF HiperGator queue system. See pca_outl_check.job for more details.
# Resources: 2Mb, 30 sec

module load plink/1.90b3.39 

INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
BASE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca"

plink --vcf "$INFILE00" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract "$BASE_DIR"/maf0.00/all_maf0.00.prune.in --vcf-half-call m --make-bed --pca --remove to_remove.fam --out all_maf0.00_outl_check
plink --vcf "$INFILE05" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract "$BASE_DIR"/maf0.05/all_maf0.05.prune.in --vcf-half-call m --make-bed --pca --remove to_remove.fam --out all_maf0.05_outl_check
```

### Plot PCAs

```R
library(ggplot2)

# MAF = 0.00
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/PCA/outlier_check"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf0.00_outl_check.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf0.00_outl_check.eigenval", sep = "/")

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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00 outlier check")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.4, 0.7)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00 outlier check, PC1 vs. PC2")
b <- b + geom_text(data = subset(pca, sample == "WG04"), aes(label = sample), hjust = -0.25, vjust = 0.25)
b

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.4, 0.7)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 outlier check, PC1 vs. PC3")
c

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00 outlier check, PC2 vs. PC3")
d <- d + geom_text(data = subset(pca, sample %in% c("WG04", "WH04")), aes(label = sample), hjust = -0.25, vjust = 0.25)
d

# MAF = 0.05
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf0.05_outl_check.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf0.05_outl_check.eigenval", sep = "/")

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
e <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
e + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.05 outlier check")
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
f <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.25, 0.375)
f <- f + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
f <- f + theme_light() + coord_equal()
f <- f + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05 outlier check, PC1 vs. PC2")
f

# PCA 1 & 3
g <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.4, 0.5)
g <- g + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
g <- g + coord_equal() + theme_light()
g <- g + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 outlier check, PC1 vs. PC3")
g

# PCA 2 + 3
h <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + xlim(-0.45, 0.4)
h <- h + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
h <- h + coord_equal() + theme_light()
h <- h + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05 outlier check, PC2 vs. PC3")
h <- h
h
```

**Results:**

After excluding the outgroup and outlier WF03/1051, plots looked fairly similar to those retaining those samples.

![Plot of PC1 vs PC2 without WF03/1051; same patterns observed as MAF=0.05 PC1 vs. PC2 with outlier sample.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/outlier_check/maf0.05_outl_pc12.png "PC1 vs. PC2")

![Plot of PC1 vs PC3; ; same patterns observed as MAF=0.05 PC1 vs. PC3 with outlier sample.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/outlier_check/maf0.05_outl_pc13.png "PC1 vs. PC3")

![Plot of PC2 vs PC3; ; same patterns observed as MAF=0.05 PC2 vs. PC3 with outlier sample.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/outlier_check/maf0.05_outl_pc23.png "PC2 vs. PC3")
