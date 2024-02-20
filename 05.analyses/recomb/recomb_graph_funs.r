# R functions to graph results of ReLERNN by chromosome and window
get_bounds <- function(tab){
    y_upper <- max(tab$CI95HI) + (0.10 * max(tab$CI95HI))
    y_lower <- min(tab$CI95LO) - (0.10 * max(tab$CI95HI))
    x_upper <- max(tab$end)
    x_lower <- 0
    return(list(xl = x_lower, xu = x_upper, yl = y_lower, yu = y_upper))
}

get_avgs <- function(tab){
    table_winsize <- unlist(apply(tab, 1, function(x) as.numeric(x["end"]) - as.numeric(x["start"])))
    avg_winsize <- mean(table_winsize)
    avg_nsites <- mean(tab$nSites)
    avg_recomb <- mean(tab$recombRate)
    return(list(winsize = avg_winsize, nsites = avg_nsites, recomb = avg_recomb))
}

graph_rwin <- function(tab, outname){
    # get plot attributes
    bounds <- get_bounds(tab)
    xlim_vec <- c(bounds[["xl"]], bounds[["xu"]])
    ylim_vec <- c(bounds[["yl"]], bounds[["yu"]])
    plot_title <- paste(tab[1, "chrom"], "Recombination Rate")
    xlab_str <- "position (bp)"
    ylab_str <- "recomb rate (per gen bp)"

    # make plot
    png(outname, width = 1400, height = 1200)

    # make plotting space
    plot(1, type = "n", main = plot_title, xlab = xlab_str, ylab = ylab_str, xlim = xlim_vec, ylim = ylim_vec, cex.main = 4, cex.lab = 2)
    
    # get midpoint of windows as x values for plotting
    recomb_pos <- unlist(apply(tab, 1, function(x) mean(c(as.numeric(x["start"]), as.numeric(x["end"])))))

    # plot conf intervals
    for(i in c(1:nrow(tab))){
        # plot conf interval
        x_points <- c(recomb_pos[i], recomb_pos[i])
        y_points <- c(tab[i, "CI95LO"], tab[i, "CI95HI"])
        lines(x = x_points, y = y_points, col = "gray", lwd = 2)
    }
    # plot calculated recombination rates
    points(x = recomb_pos, y = tab$recombRate, pch = 19, col = "black")

    # plot chromosome-wide averages
    avgs <- get_avgs(tab)
    avg_txt <- paste("avg=", as.character(signif(avgs[["recomb"]], 3), "/gen*bp", sep = ""))
    # get coordinates for top right corner of plot
    tr_x <- bounds[["xu"]] - (0.02 * bounds[["xu"]])
    tr_y <- bounds[["yu"]] - (0.02 * bounds[["yu"]])
    text(tr_x, tr_y, labels = avg_txt, adj = c(1, 1), cex = 4)

    dev.off()
}