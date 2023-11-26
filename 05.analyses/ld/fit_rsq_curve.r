# R script to fit decay curve to SNP correlation values, to derive linkage disequilibrium from emeraLD output
# Usage: Rscript fit_rsq_curve.r [rsq_table] [output_name] [min_comps] [r2_cutoff] [graph_bool]
#        rsq_table: average r2 value per distance between variants from emeraLD with the columns dist,r2
#        output_name: name of output table to print
#        min_comps: minimum number of comparisons made in a window to include it in output
#        r2_cutoff: r2 at which to report distance for linkage disequilibrium measure
#        graph_bool: boolean to turn on (1) or off (0) the printing of graphs of r2 vs. distance

# Output is a table with the filename, distance at r2 = cutoff, estimated C from model
# Uses LD decay model from Hill and Weiss 1988

# Import libraries, get arguments
file_list_name <- commandArgs(trailingOnly = TRUE)[1]
out_name <- commandArgs(trailingOnly = TRUE)[2]
min_comps <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
r2_cutoff <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
graph_bool <- as.logical(commandArgs(trailingOnly = TRUE)[5])

# Define model
exp_rsquare <- function(l, C) (((10+C*l)/((2+C*l)*(11+C*l)))*(1+(((3+C*l)*(12+12*C*l+(C*l)^2))/(30*(2+C*l)*(11+C*l)))))
# simple_rsquare <- function(l, C) 1/(1+C*l)

# Import r-square tables
tab_list <- read.table(file_list_name, header = FALSE)$V1
# Initialize output table
out_tab <- data.frame(file_name = character(), ld = numeric(), C = numeric())

for(filename in tab_list){
    stripped_filename <- unlist(strsplit(filename, "/"))[length(unlist(strsplit(filename, "/")))]
    # Import file
    print(paste("doing", filename))
    rsq_tab <- read.csv(filename, header = TRUE)
    rsq_tab <- rsq_tab[which(rsq_tab$r2 != -1),]
    # Check that it has enough observations to estimate LD
    if(nrow(rsq_tab) > min_comps){
        # Fit model
        rsq_model <- nls(r2 ~ exp_rsquare(dist, estC), data = rsq_tab, start = list(estC = 1))
        print(summary(rsq_model))
        # calculate approximate LD at r2 cutoff by interpolated from predicted curve
        est_C <- rsq_model$m$getPars()["estC"]
        pred_dists <- c(1:max(rsq_tab$dist))
        pred_curve <- exp_rsquare(pred_dists, est_C)
        above_cutoff <- max(which(pred_curve > r2_cutoff)) # returns a bp
        below_cutoff <- min(which(pred_curve < r2_cutoff)) # returns a bp
        cutoff_ld <- mean(above_cutoff, below_cutoff)
        # store estimated C and LD at r2 cutoff
        out_tab <- rbind(out_tab, list(file_name = stripped_filename, ld = cutoff_ld, C = est_C))
        # Only plot if switch is turned on
        if(graph_bool){
            # Plot model against actual data points
            png(paste(stripped_filename, "_curve.png", sep = ""))
            plot(rsq_tab$dist, rsq_tab$r2, main = paste(stripped_filename, "MAC = 1 LD"), xlab = "distance (bp)", ylab = "r2")
            lines(rsq_tab$dist, predict(rsq_model), col = "cyan")
            dev.off()
            # Plot residuals
            png(paste(stripped_filename, "_resid.png", sep = ""))
            plot(rsq_tab$dist, residuals(rsq_model), main = paste(stripped_filename, "MAC = 1 LD residuals"))
            dev.off()
        }
    } else {
        print(paste(filename, "did not have enough observations"))
    }
    # Export output table
    write.table(out_tab, out_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}