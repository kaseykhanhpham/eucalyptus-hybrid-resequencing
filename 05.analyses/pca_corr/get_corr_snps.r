# R script to retrieve correlated SNPs to a PC using PCAdapt
# Usage: Rscript get_corr_snps.r [bedfile] [bimfile] [k] [alpha]

library(pcadapt)
library(qvalue)

bedfile_name <- commandArgs(trailingOnly = TRUE)[1]
bimfile_name <- commandArgs(trailingOnly = TRUE)[2]
k <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
alpha <- as.numeric(commandArgs(trailingOnly = TRUE)[4])

# get correlated SNPs to first k PCs using main pcadapt function
bed_in <- read.pcadapt(input = bedfile_name, type = "bed")
bim_in <- read.table(bimfile_name, sep = "\t", header = FALSE)
corr_test <- pcadapt(input = bed_in, K = k)

# print results of pcadapt correlation test
write(summary(corr_test), stdout())

# export diagnostic plots
file_only <- unlist(strsplit(bedfile_name, "/"))[length(unlist(strsplit(bedfile_name, "/")))]
name_root <- gsub(".bed$", "", file_only)

png(paste(name_root, "_manhattan.png", sep = ""), width = 1000, height = 1000, units = "px")
    plot(corr_test, option = "manhattan")
dev.off()
png(paste(name_root, "_qqplot.png", sep = ""), width = 1000, height = 1000, units = "px")
    plot(corr_test, option = "qqplot")
dev.off()
png(paste(name_root, "_pvals.png", sep = ""), width = 1000, height = 1000, units = "px")
    hist(corr_test$pvalues, xlab = "p-values", main = NULL, breaks = 50)
dev.off()
png(paste(name_root, "_statdistr.png", sep = ""), width = 1000, height = 1000, units = "px")
    plot(corr_test, option = "stat.distribution")
dev.off()

## Get outliers
qval <- qvalue(corr_test$pvalues)$qvalues 
outliers <- which(qval < alpha)

bim_out <- bim_in[outliers,]
write.table(bim_out, paste(name_root, "_outliers.bim", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")