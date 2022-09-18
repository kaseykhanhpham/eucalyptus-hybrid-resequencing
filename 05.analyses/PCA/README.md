# Whole-genome PCA
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/pca/).

### Prune linked SNPs
Done for both the maf=0.00 and maf=0.05 SNP sets since I use both extensively in other analyses. Note: this step will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see link_prune.job for details
# Resources used: 1 Gb, 30 sec

module load plink/1.90b3.39 

export INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00_snps.vcf"
export INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05_snps.vcf"

plink --vcf "$INFILE00" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --vcf-half-call m --out all_maf0.00
plink --vcf "$INFILE05" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --vcf-half-call m --out all_maf0.05

```

### Calculate principal components
Note: this step also will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see pca.job for details
# Resources used: 2 Mb, 9 sec

module load plink/1.90b3.39 

export INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05_snps.vcf"

plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract all_maf0.05.prune.in --vcf-half-call m --make-bed --pca --out all_maf0.05
```

### Plot principal components
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

Two samples, WF03 and WG04, were very different from the rest of their groups (reference _E. globulus_ and introgressed _E. globulus_ respectively) for MAF=0.00 plots. WG04 fell back within range of the rest of introgressed _E. globulus_ variation in MAF=0.05 plots. Below, MAF=0.05 plots are shown.


![Plot of PC1 vs PC2; PC1 differentiates strongly between _E. cordata_ and _E. globulus_ (both introgressant and pure), while PC2 differentiates between introgressant and pure _E. globulus_, with _E. cordata_ clustering with introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc12_labeled.png "PC1 vs. PC2")

![Plot of PC1 vs PC3; very similar to PC1 vs PC2 but _E. globulus_ populations intergrade more along PC3.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc13_labeled.png "PC1 vs. PC3")

![Plot of PC2 vs PC3; pure _E. globulus_ seems to cluster slightly away from the other individuals, while _E. cordata_ clusters tightly in the center of a cloud of introgressant _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc23_labeled.png "PC2 vs. PC3")

![Barplot of percentage genetic variation explained by each principal component; PC1 explains about 30% of variation, while all other PCs pictured explain about 5%](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/PCA/maf0.05/maf05_pc_var_explained.png "Percent Variance Explained by each PC")