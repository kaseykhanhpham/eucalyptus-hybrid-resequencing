prefix <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(dplyr)

print(paste("doing ", prefix, sep = ""))

# log(QUAL)
var_qual <- read.table(paste(prefix, ".lqual", sep = ""), header = TRUE)
logqual <- log(var_qual$QUAL)
var_qual <- cbind(var_qual, logqual)
a <- ggplot(var_qual, aes(logqual)) + geom_density(fill = "grey", colour = "black")
png(paste(prefix, "_logqual.png", sep = ""))
a + theme_light()
dev.off()

# DEPTH
var_depth <- read.table(paste(prefix, ".ldepth.mean", sep = ""), header = TRUE)
a <- ggplot(var_depth, aes(MEAN_DEPTH)) + geom_density(fill = "grey", colour = "black")
png(paste(prefix, "_dpbysite.png", sep = ""))
a + theme_light()
dev.off()

# MISSING
var_miss <- read.table(paste(prefix, ".lmiss", sep = ""), header = TRUE)
a <- ggplot(var_miss, aes(F_MISS)) + geom_density(fill = "grey", colour = "black")
png(paste(prefix, "_missing.png", sep = ""))
a + theme_light()
dev.off()

# IND MISSING
ind_miss  <- read.table(paste(prefix, ".imiss", sep = ""), header = TRUE)
a <- ggplot(ind_miss, aes(F_MISS)) + geom_histogram(fill = "grey", colour = "black")
png(paste(prefix, "_samplemissing.png", sep = ""))
a + theme_light()
dev.off()

# MAF
var_freq <- read.table(paste(prefix, ".frq", sep = ""), header = TRUE)
MAF <- apply(var_freq[, c("A1", "A2")], 1, function(z) min(z))
var_freq <- cbind(var_freq, MAF)
a <- ggplot(var_freq, aes(MAF)) + geom_density(fill = "grey", colour = "black")
png(paste(prefix, "_maf.png", sep = ""))
a + theme_light()
dev.off()

# HETEROZYGOSITY
ind_het <- read.table(paste(prefix, ".het", sep = ""), header = TRUE)
a <- ggplot(ind_het, aes(F)) + geom_histogram(fill = "grey", colour = "black")
png(paste(prefix, "_het.png", sep = ""))
a + theme_light()
dev.off()

