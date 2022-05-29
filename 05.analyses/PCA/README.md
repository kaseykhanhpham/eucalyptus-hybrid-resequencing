# Whole-genome PCA
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/pca/).

### Prune linked SNPs
Note: this step will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see link_prune.job for details
# Resources used: 2 Mb, 13 sec

module load plink/1.90b3.39 
export INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05.vcf"

plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --vcf-half-call m --out all_maf0.05
```

### Calculate principal components
Note: this step also will not work as given in my job files if using `plink v2` or higher.

```bash
# Run via job on UFRC, see pca.job for details
# Resources used: 2 Mb, 9 sec

module load plink/1.90b3.39 

export INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05.vcf"

plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract all_maf0.05.prune.in --vcf-half-call m --make-bed --pca --out all_maf0.05
```

### Plot principal components
Performed locally in [`R programming language`](https://www.r-project.org/).

```R
library(ggplot2)
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/euc_hyb_reseq/Analyses/PCA/"
table_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/euc_hyb_reseq/sample_spp_table.csv"
pca_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/euc_hyb_reseq/Analyses/PCA/all_maf0.05.eigenvec"
eigenval_address <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/euc_hyb_reseq/Analyses/PCA/all_maf0.05.eigenval"

# Read input files
setwd(working_dir)
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
a + ylab("Proportion Variance Explained") + theme_light()
cumsum(pve$pve)

# Plot PCA
# PCA 1 & 2
b <- ggplot(pca, aes(PC1, PC2, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.6)
b <- b + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
b <- b + theme_light() + coord_equal()
b <- b + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = ""))
b
# PCA 1 & 3
c <- ggplot(pca, aes(PC1, PC3, col = taxon)) + geom_point(size = 2) + xlim(-0.2, 0.8)
c <- c + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
c <- c + coord_equal() + theme_light()
c <- c + xlab(paste("PC1 (", signif(pve$pve[1], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = ""))
c
# PCA 2 + 3
d <- ggplot(pca, aes(PC2, PC3, col = taxon, label = sample)) + geom_point(size = 2) + xlim(-0.30, 0.65) + geom_text(hjust=-.1,vjust=.1)
d <- d + scale_colour_manual(values = c("goldenrod1", "black", "deepskyblue4"))
d <- d + coord_equal() + theme_light()
d <- d + xlab(paste("PC2 (", signif(pve$pve[2], 3), ")", sep = "")) + ylab(paste("PC3 (", signif(pve$pve[3], 3), ")", sep = ""))
d
```