prefix <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(dplyr)

print(paste("doing ", prefix, sep = ""))

# log(QUAL)
var_qual <- read.table(paste(prefix, ".lqual", sep = ""), header = TRUE)
logqual <- log(var_qual$QUAL)
print("QUALITY")
print(summary(var_qual$QUAL))
print("LOGQUALITY")
print(summary(logqual))

# DEPTH
var_depth <- read.table(paste(prefix, ".ldepth.mean", sep = ""), header = TRUE)
print("DEPTH")
print(summary(var_depth$MEAN_DEPTH))

# MISSING
var_miss <- read.table(paste(prefix, ".lmiss", sep = ""), header = TRUE)
print("MISSING")
print(summary(var_miss$F_MISS))

# MAF
var_freq <- read.table(paste(prefix, ".frq", sep = ""), header = TRUE)
MAF <- apply(var_freq[, c("A1", "A2")], 1, function(z) min(z))
print("MAF")
print(summary(MAF))

