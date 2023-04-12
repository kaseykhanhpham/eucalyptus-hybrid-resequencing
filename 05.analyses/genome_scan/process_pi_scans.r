# R script to process egglib pi genome scan output
# Usage: Rscript process_pi_scans.r [infile]

infile_name <- commandArgs(trailingOnly=TRUE)[1]
pi_table <- read.table(infile_name, header = TRUE, sep = "\t", na.strings="")

print(paste("Sample", infile_name, ":"))
print(summary(pi_table$Pi, na.rm = TRUE))
png(paste(infile_name, ".pi_distr.png", sep = ""), width = 500, height = 500)
    hist(pi_table$Pi)
dev.off()

pi_mean <- mean(pi_table$Pi, na.rm = TRUE)
pi_sd <- sd(pi_table$Pi, na.rm = TRUE)

flagged_pi <- which(pi_table$Pi > (mean(pi_mean + 2*pi_sd)))
flagged_pi_windows <- pi_table[flagged_pi, c("chr", "start", "end")]

write.table(flagged_pi_windows, paste(infile_name, ".flagged_pi.tab", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE)