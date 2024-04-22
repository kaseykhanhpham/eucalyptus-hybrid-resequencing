# Whole-genome PCA
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/pca/).

## Prune linked SNPs
Note: this step will not work as given in my job files if using `plink v2` or higher.

```bash
# Run on UFRC queue system; see link_prune.job for details.
# Resources used: 15 Mb, 2 min

module load plink/1.90b3.39 

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
OUTDIR00="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00"
OUTDIR05="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    plink --vcf "$INDIR"/"$NAME"_fil.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --vcf-half-call m --out "$OUTDIR00"/"$NAME"_maf00
    plink --vcf "$INDIR"/"$NAME"_fil.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --vcf-half-call m --maf 0.05 --out "$OUTDIR05"/"$NAME"_maf05
done
```

Merge pruned chromosome PLINK files.

```bash
OUTDIR00="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf00"
OUTDIR05="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf05"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

touch "$OUTDIR00"/all_maf00.prune.in
touch "$OUTDIR05"/all_maf05.prune.in

for NAME in "${VCFLIST[@]}"
do
    cat "$OUTDIR00"/"$NAME"_maf00.prune.in >> "$OUTDIR00"/all_maf00.prune.in
    cat "$OUTDIR05"/"$NAME"_maf05.prune.in >> "$OUTDIR05"/all_maf05.prune.in
done
```

Calculate PCAs in `plink`.

```bash
# Run in UFRC queue system; see pca.job for more details.
# Resources used: 1 Gb, 4 min

module load plink/1.90b3.39 

VCF="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"
OUTDIR00="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf00"
OUTDIR05="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf05"

plink --vcf "$VCF" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract "$OUTDIR00"/all_maf00.prune.in --vcf-half-call m --make-bed --pca --out "$OUTDIR00"/all_maf00
plink --vcf "$VCF" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract "$OUTDIR05"/all_maf05.prune.in --vcf-half-call m --maf 0.05 --make-bed --pca --out "$OUTDIR05"/all_maf05
```

## Plot principal components
Performed locally in [`R programming language`](https://www.r-project.org/). Used metadata to color samples by type: reference _E. cordata_ ("cord_MR", yellow), reference _E. globulus_ ("glob_ref", blue), or introgressed _E. globulus_ ("glob_MR", black).

### MAF=0.00

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/wg_pca/maf0.00"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf00.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf00.eigenval", sep = "/")

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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.00, ingroup only")
cumsum(pve$pve)

# Scale by percent variation
pca_scaled <- as.data.frame(sapply(c(1:20), function(i) pve[i, "pve"] * pca[, paste("PC", i, sep = "")]))
colnames(pca_scaled) <- colnames(pca)[2:21]
pca_scaled$taxon <- pca$taxon

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca_scaled, aes(PC1, PC2, col = taxon)) + geom_point(size = 2)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.00, PC1 vs. PC2")
b <- b + ylim(-0.075, 0.075)
b

# PCA 1 & 3
c <- ggplot(pca_scaled, aes(PC1, PC3, col = taxon)) + geom_point(size = 2)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00, PC1 vs. PC3")
c <- c + ylim(-0.075, 0.075) 
c

# PCA 2 + 3
d <- ggplot(pca_scaled, aes(PC2, PC3, col = taxon)) + geom_point(size = 2)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.00, PC2 vs. PC3")
d <- d + ylim(-0.02, 0.02)
d
```

### MAF=0.05

```R
library(ggplot2)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/wg_pca/maf0.05"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
pca_address <- paste(working_dir, "all_maf05.eigenvec", sep = "/")
eigenval_address <- paste(working_dir, "all_maf05.eigenval", sep = "/")

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
a + ylab("Proportion Variance Explained") + theme_light() + ggtitle("maf=0.05 Eigenvalues")
cumsum(pve$pve)

# Scale PCs to eigenvalue
pca_scaled <- as.data.frame(sapply(c(1:20), function(i) pve[i, "pve"] * pca[, paste("PC", i, sep = "")]))
colnames(pca_scaled) <- colnames(pca)[2:21]
pca_scaled$taxon <- pca$taxon

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ggtitle("maf=0.05, PC1 vs. PC2")
b <- b + xlim(-0.4, 0.4)
b

# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05, PC1 vs. PC3")
c <- c + xlim(-0.45, 0.50)
c

# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon)) + geom_point(size = 2)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = "")) + ggtitle("maf=0.05, PC2 vs. PC3")
d <- d + xlim(-0.45, 0.5)
d
```

## Results


**Overall patterns:** PC1 explained the vast majority of variation (~30%) for both the MAF=0.00 and MAF=0.05 datasets, with all other PCs explaining ~5% of present variation. In both datasets and regardless of whether the outgroup was included, PC1 strongly differentiated between _E. cordata_ and reference and Meehan Range _E. globulus_. For MAF=0.00, _E. globulus_ samples continued to group together along PC2, but _E. cordata_ was split into two variable clusters (likely the family clusters). This was not apparent in the MAF=0.05 dataset, where _E. cordata_ clustered together tightly between the _E. globulus_ reference group and the Meehan Range _E. globulus_ reference group. This indicates that for these filtering parameters, _E. cordata_ samples contain a lot of singleton variation separating the two families.

### MAF = 0.00 Plots

![Plot of PC1 vs PC2 for MAF=0.00 (excluding the outgroup); PC1 differentiates strongly between _E. cordata_ and _E. globulus_ (both introgressant and pure), while _E. cordata_ is split into two groups along PC2.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.00/maf00_pc12.png "MAF=0.00 PC1 vs. PC2")

![Plot of PC1 vs PC3 for MAF=0.00 (excluding the outgroup); _E. globulus_ samples differentiate along PC3 with _E. cordata_ sitting between the two, closer to Meehan Range _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.00/maf00_pc13.png "MAF=0.00 PC1 vs. PC3")

![Plot of PC2 vs PC3 for MAF=0.00 (excluding the outgroup)](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.00/maf00_pc23.png "MAF=0.00 PC2 vs. PC3")

![Barplot of percentage genetic variation explained by each principal component; PC1 explains about 30% of variation, while all other PCs pictured explain about 5%](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.00/maf00_var_expl.png "MAF=0.00 Percent Variance Explained by each PC")

### MAF = 0.05 Plots

![Plot of PC1 vs PC2 for MAF=0.05 (exlcluding the outgroup); PC1 differentiates strongly between _E. cordata_ and _E. globulus_ (both introgressant and pure), while _E. cordata_ associates slightly more with Meehan Range _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.05/maf05_pc12.png "MAF=0.05 PC1 vs. PC2")

![Plot of PC1 vs PC3 for MAF=0.05 (excluding the outgroup); Sample WF03/1051 separates from the other samples along PC3.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.05/maf05_pc13.png "MAF=0.05 PC1 vs. PC3")

![Plot of PC2 vs PC3 for MAF=0.00 (excluding the outgroup)](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.05/maf05_pc23.png "MAF=0.05 PC2 vs. PC3")

![Barplot of percentage genetic variation explained by each principal component; PC1 explains about 30% of variation, while all other PCs pictured explain about 5%](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_pca/maf0.05/maf05_var_expl.png "MAF=0.05 Percent Variance Explained by each PC")
