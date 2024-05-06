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
    # dependent on ggplot2 package

    bounds <- get_bounds(tab)
    # xlim_vec <- c(bounds[["xl"]], bounds[["xu"]])
    # ylim_vec <- c(bounds[["yl"]], bounds[["yu"]])
    plot_title <- paste(tab[1, "chrom"], "Recombination Rate")
    xlab_str <- "position (bp)"
    ylab_str <- "recomb rate (per gen bp)"

    # get midpoint of windows as x values for plotting
    recomb_pos <- unlist(apply(tab, 1, function(x) mean(c(as.numeric(x["start"]), as.numeric(x["end"])))))
    tab$recomb_pos <- recomb_pos

    # make plot
    png(outname, width = 1400, height = 1200)

    p <- ggplot(data = tab, aes(x = recomb_pos, y = recombRate)) + geom_point(size = 4) + geom_line() + theme_bw(base_size = 28) + xlab(xlab_str) + ylab(ylab_str) + ggtitle(plot_title)
    p <- p + geom_ribbon(aes(ymin = CI95LO, ymax = CI95HI), linetype = 2, alpha = 0.35)

    # plot chromosome-wide averages
    avgs <- get_avgs(tab)
    avg_txt <- paste("avg =", as.character(signif(avgs[["recomb"]], 3), "/gen*bp", sep = ""))
    # get coordinates for top right corner of plot
    tr_x <- bounds[["xu"]] - (0.02 * bounds[["xu"]])
    tr_y <- bounds[["yu"]] - (0.02 * bounds[["yu"]])
    p <- p + annotate("text", x = tr_x, y = tr_y, label = avg_txt, adj = c(1, 1), size = 12)

    print(p)

    dev.off()
}