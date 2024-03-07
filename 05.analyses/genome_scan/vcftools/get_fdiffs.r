# Rscript to extract fixed differences from output of vcftools site-pi for
# E. globulus (reference) and E. cordata

# Import files
glob_tab <- read.table("globref.sites.pi", header = TRUE, sep = "\t")
cord_tab <- read.table("cord.sites.pi", header = TRUE, sep = "\t")
refs_tab <- read.table("refs.sites.pi", header = TRUE, sep = "\t")

# Create storage variables for loops
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11", "ChrUn")
chr_full <- c()
sites_full <- c()
for(chr in chr_list){
    # subset tables to just the current chromosome
    glob_fil <- glob_tab[which(glob_tab$CHROM == chr),]
    cord_fil <- cord_tab[which(cord_tab$CHROM == chr),]
    refs_fil <- refs_tab[which(refs_tab$CHROM == chr),]
    # Get full list of sites on the chromosome that appear in any of the tables
    chr_sites <- order(unique(c(glob_fil$POS, cord_fil$POS, refs_fil$POS)))
    sites_include <- c()
    for(site in chr_sites){
        # Check for presence of site in each table and if present, the value
        # looking for 0 in glob and cord only and != 0 in both table
        if(site %in% glob_fil$POS){
            glob_status <- glob_fil[which(glob_fil$POS == site), "PI"] == 0
        } else {
            glob_status <- FALSE
        }
        if(site %in% cord_fil$POS){
            cord_status <- cord_fil[which(cord_fil$POS == site), "PI"] == 0
        } else {
            cord_status <- FALSE
        }
        if(site %in% refs_fil$POS){
            refs_status <- refs_fil[which(refs_fil$POS == site), "PI"] != 0
        } else {
            refs_status <- FALSE
        }
        # add site to final list if 0 in glob and cord tables and != 0 in table with both
        if(glob_status & cord_status & refs_status){
            sites_include <- c(sites_include, site)
        }
    }
    # after looping through all sites on a chromosome, append the list to the genome-wide list
    chr_full <- c(chr_full, rep(chr, length(sites_include)))
    sites_full <- c(sites_full, sites_include)
}
# Write final list of fixed differences to output file
fdiffs_tab <- data.frame(chr = chr_full, pos = sites_full)
write.table(fdiffs_tab, "fixed_diffs_globref_cord.tab", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")