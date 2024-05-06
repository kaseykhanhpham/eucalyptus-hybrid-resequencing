# Synthesis analysis and plotting
## Identifying Regions of Interest
Retrieved regions of interest and plotted chromosome painting-style by synthesizing the genome scan stats, sliding window LD, sliding window recombination, and local ancestry inference results. Refer to `win_processing_funs.r`, `plot_chr_wins.r`, `intr_lais_funs.r`, `plot_intr_lais.r` for details on process.

Regions of interest for species differences, recombination, introgression, and selection were defined as follows:

### Species Differences

**Populations:** Between pure _E. globulus_ and _E. cordata_.
**Minimum sites to consider a window:** 15

| Category        | FST cutoff      | site pi         | 
| --------------- | --------------- | --------------- |
| Moderate Diffs  | 80th percentile | None            |
| Large Diffs     | 95th percentile | None            |
| Fixed Diffs     | None            | pi = 0 within species and pi > 0 between species |

### Recombination

**Populations:** Meehan Range _E. globulus_

| Category        | Recombination cutoff | LD cutoff       |
| --------------- | -------------------- | --------------- |
| Hotspots        | 95th percentile      | None            |
| Suppressed      | Lowest 5% of values  | 80th percentile |

### Introgression

**Populations:** Between Meehan Range _E. globulus_ and _E. cordata_.
**Minimum sites to consider a window:** 40
**Minimum individuals sharing LAI score to consider a window:** 5

| Category        | dxy cutoff           | fdM cutoff      | LAI cutoff        |
| --------------- | -------------------- | --------------- | ----------------- |
| Genome scan     | Lowest 40% of values | 90th percentile | None              |
| AHMM regions    | None                 | None            | Posterior > 0.95  |
| ELAI regions    | None                 | None            | Avg Dosage > 1.75 |

### Selection

**Populations:** Meehan Range _E. globulus_
**Minimum sites to consider a window:** 40

| Category              | Tajima's D cutoff   | Recombination cutoff | LD cutoff        |
| --------------------- | ------------------- | -------------------- | ---------------- |
| Balancing selection   | 95th percentile     | 50th percentile      | 75th percentile  |
| Directional selection | Lowest 5% of values | 50th percentile      | 75th percentile  |

## Synthesis
Wanted to know if species difference loci overlap more with selection loci than expected by chance. First converted FST output from `pixy` into a BED file.

```R
fst_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/all_fst.txt"
fst_file <- read.table(fst_filename, header = TRUE, sep = "\t")
fst_MRcord_bed <- fst_file[which(fst_file$pop1 == "glob_MR" & fst_file$pop2 == "cord_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
colnames(fst_MRcord_bed) <- c("chrom", "start", "end")
fst_purecord_bed <- fst_file[which(fst_file$pop1 == "glob_pure" & fst_file$pop2 == "cord_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
colnames(fst_purecord_bed) <- c("chrom", "start", "end")
write.table(fst_MRcord_bed, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/globMR_cordMR_fst.bed", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(fst_purecord_bed, "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/globpure_cordMR_fst.bed", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

Got intersection of the two categories in `BedTools` (strong species differences category).

```bash
module load bedtools/2.30.0
FDIFFS_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/fdiff_files"
SEL_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files"
FST_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"

# balancing
bedtools intersect -a "$FDIFFS_LOC"/Eglob_all_fdiffs.bed -b "$SEL_LOC"/balancing/Eglob_all_balance_sel.bed -header > Eglob_fdiffs_balsel_inters.bed
# directional
bedtools intersect -a "$FDIFFS_LOC"/Eglob_all_fdiffs.bed -b "$SEL_LOC"/directional/Eglob_all_direction_sel.bed -header > Eglob_fdiffs_dirsel_inters.bed
# merge balancing and directional BED files
cat "$SEL_LOC"/balancing/Eglob_all_balance_sel.bed > "$SEL_LOC"/Eglob_all_both_sel.bed
cat "$SEL_LOC"/directional/Eglob_all_direction_sel.bed >> "$SEL_LOC"/Eglob_all_both_sel.bed
# get FST windows between E. glob (MR) and E. cord at least 75% overlapping with both types of selection windows
bedtools intersect -a "$FST_LOC"/globMR_cordMR_fst.bed -b "$SEL_LOC"/Eglob_all_both_sel.bed -u -f 0.75 -header > globMR_cordMR_fst_sel.bed
# get FST windows between E. glob (pure) and E. cord at least 75% overlapping with selection windows
bedtools intersect -a "$FST_LOC"/globpure_cordMR_fst.bed -b "$SEL_LOC"/Eglob_all_both_sel.bed -u -f 0.75 -header > globpure_cordMR_fst_sel.bed
```

Calculated coverage across genome of each set of regions separately and their overlap in `R v4.2`.
```R
## Compare overlap between species difference and selection windows ##
#  ----------------------------------------------------------------- #
# Import BED files as data frames
fdiffs_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/fdiff_files/Eglob_all_fdiffs.bed"
balsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/balancing/Eglob_all_balance_sel.bed"
dirsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/directional/Eglob_all_direction_sel.bed"
balsel_overl_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/Eglob_fdiffs_balsel_inters.bed"
dirsel_overl_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/Eglob_fdiffs_dirsel_inters.bed"

fdiffs_file <- read.table(fdiffs_filename, header = TRUE, sep = "\t")
balsel_file <- read.table(balsel_filename, header = TRUE, sep = "\t")
dirsel_file <- read.table(dirsel_filename, header = TRUE, sep = "\t")
balsel_overl_file <- read.table(balsel_overl_filename, header = TRUE, sep = "\t")
dirsel_overl_file <- read.table(dirsel_overl_filename, header = TRUE, sep = "\t")

genome_size <- 603301446

# Get total base pairs encompassed over the genome for each category
fdiffs_bp <- sum(unlist(apply(fdiffs_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
balsel_bp <- sum(unlist(apply(balsel_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
dirsel_bp <- sum(unlist(apply(dirsel_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
balsel_overl_bp <- sum(unlist(apply(balsel_overl_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
dirsel_overl_bp <- sum(unlist(apply(dirsel_overl_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))

# Get fraction of whole genome for each category
fdiffs_frac <- fdiffs_bp/genome_size
balsel_frac <- balsel_bp/genome_size
dirsel_frac <- dirsel_bp/genome_size
balsel_overl_frac <- balsel_overl_bp/genome_size
dirsel_overl_frac <- dirsel_overl_bp/genome_size

# Get expected balancing selection overlap and directional selection overlap expected values
balsel_overl_expected <- fdiffs_frac * balsel_frac * genome_size
dirsel_overl_expected <- fdiffs_frac * dirsel_frac * genome_size

## Get average FST of windows identified as under balancing selection or directional selection ##
# --------------------------------------------------------------------------------------------- #
# pixy output with window fst values
fst_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/all_fst.txt"
fst_file <- read.table(fst_filename, header = TRUE, sep = "\t")
# BED files of just windows for each population pair
fst_sel_MRcord_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/globMR_cordMR_fst_sel.bed"
fst_sel_MRcord_bed <- read.table(fst_sel_MRcord_filename, header = TRUE, sep = "\t")
fst_sel_purecord_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/globpure_cordMR_fst_sel.bed"
fst_sel_purecord_bed <- read.table(fst_sel_purecord_filename, header = TRUE, sep = "\t")

min_snps <- 15 # set minimum SNP number to consider a window

# define filters
# snp number filter
snp_fil <- which(fst_file$no_snps > min_snps)
# filter for correct population pairs
fst_MRcord_popfil <- which(fst_file$pop1 == "glob_MR" & fst_file$pop2 == "cord_MR")
fst_purecord_popfil <- which(fst_file$pop1 == "glob_pure" & fst_file$pop2 == "cord_MR")
# get rows that are also overlapping with selection windows, as defined in imported BED file
sel_MRcord_winfil <- unlist(apply(fst_sel_MRcord_bed, 1, function(row) which(fst_file$chromosome == row["chrom"] & fst_file$window_pos_1 == as.integer(row["start"]) & fst_file$window_pos_2 == as.integer(row["end"]))))
sel_purecord_winfil <- unlist(apply(fst_sel_purecord_bed, 1, function(row) which(fst_file$chromosome == row["chrom"] & fst_file$window_pos_1 == as.integer(row["start"]) & fst_file$window_pos_2 == as.integer(row["end"]))))
# final filters
MRcord_final_fil <- intersect(snp_fil, intersect(fst_MRcord_popfil, sel_MRcord_winfil))
purecord_final_fil <- intersect(snp_fil, intersect(fst_purecord_popfil, sel_purecord_winfil))

# get final filtered tables for each population pair
MRcord_fst_sel_tab <- fst_file[MRcord_final_fil,]
purecord_fst_sel_tab <- fst_file[purecord_final_fil,]

# get average FST
MRcord_fst_sel_avg <- mean(MRcord_fst_sel_tab$avg_wc_fst) # 0.3658838
purecord_fst_sel_avg <- mean(purecord_fst_sel_tab$avg_wc_fst) # 0.3471521
```

| Category                          | bp covered      | fraction of genome | expected bp covered | difference (bp)|
| --------------------------------- | --------------- | ------------------ | ------------------- | -------------- |
| species differences               | 5958976         | 0.009877278        | NA                  | NA             |
| balancing selection               | 2501955         | 0.004147106        | NA                  | NA             |
| directional selection             | 3042942         | 0.005043817        | NA                  | NA             |
| spp diffs & balancing selection   | 14997           | 0.00002485822      | 24712.5             | -9715.5        |
| spp diffs & directional selection | 197973          | 0.0003281494       | 30055.98            | 167917.02      |

The amount of balancing selection windows overlapping with species difference windows is about what you'd expect from random chance. However, the amount of directional selection windows overlapping with species difference windows is much higher, by a factor of ~6. This indicates that divergent selection on _E. globulus_ may be acting against any introgression from _E. cordata_.

As for average FST of selection windows...
| Category            | FST           | Genome-wide FST |
| ------------------- | ------------- | --------------- |
| glob_MR & cord_MR   | 0.3658838     | 0.281874        |
| glob_pure & cord_MR | 0.3471521     | 0.2748824       |

Average FST in both types of selection windows combined is also quite a bit higher than genome-wide FST, filtering for a minimum number of 15 SNPs per window (which of course biases the numbers towards regions of the genome with more coverage).

### Error Estimation
Bootstrapped selection windows to get errors on point estimates of species difference/selection window overlap and average FST in selection windows. 

#### Window Overlap
Started by generating the bootstraps in `R`.
```R
# Import BED files as data frames
balsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/balancing/Eglob_all_balance_sel.bed"
dirsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/directional/Eglob_all_direction_sel.bed"

balsel_file <- read.table(balsel_filename, header = TRUE, sep = "\t")
dirsel_file <- read.table(dirsel_filename, header = TRUE, sep = "\t")

# Resample selection windows (100x) and write resampled BED files to directories
for(i in c(1:100)){
    # generate and sort resampled indices
    balsel_res_ind <- sort(sample(c(1:nrow(balsel_file)), nrow(balsel_file), replace = TRUE))
    dirsel_res_ind <- sort(sample(c(1:nrow(dirsel_file)), nrow(dirsel_file), replace = TRUE))
    # subset original files with resampled indices
    balsel_res <- balsel_file[balsel_res_ind,] 
    dirsel_res <- dirsel_file[dirsel_res_ind,]
    # write resampled files
    write.table(balsel_res, paste("boots/Eglob_all_balance_sel/boot_", i, ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    write.table(dirsel_res, paste("boots/Eglob_all_direction_sel/boot_", i, ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
```

Got overlap of bootstrapped selection BED files with species difference windows.
```bash
module load bedtools/2.30.0
FDIFFS_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/fdiff_files"
BBOOTS_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_balance_sel"
DBOOTS_LOC="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_direction_sel"

for i in {1..100}
do
    bedtools intersect -a "$FDIFFS_LOC"/Eglob_all_fdiffs.bed -b "$BBOOTS_LOC"/boot_"$i".bed -header > boots/Eglob_all_balance_sel_overl/boot_"$i"_overl.bed
    bedtools intersect -a "$FDIFFS_LOC"/Eglob_all_fdiffs.bed -b "$DBOOTS_LOC"/boot_"$i".bed -header > boots/Eglob_all_direction_sel_overl/boot_"$i"_overl.bed
done
```

Calculated overlap bp for each bootstrap in `R`.
```R
fdiffs_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/fdiff_files/Eglob_all_fdiffs.bed"
balsel_boots_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_balance_sel"
dirsel_boots_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_direction_sel"
balsel_boots_overl_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_balance_sel_overl"
dirsel_boots_overl_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/boots/Eglob_all_direction_sel_overl"

fdiffs_file <- read.table(fdiffs_filename, header = TRUE, sep = "\t")
fdiffs_bp <- sum(unlist(apply(fdiffs_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
genome_size <- 603301446

# Calculate bootstrap estimates for overlaps
balsel_overl_vec <- c()
dirsel_overl_vec <- c()
balsel_overl_expected_vec <- c()
dirsel_overl_expected_vec <- c()

for(i in c(1:100)){
    # import boostrap files
    balsel_file <- read.table(paste(balsel_boots_dir, "/boot_", i, ".bed", sep = ""), header = TRUE, sep = "\t")
    dirsel_file <- read.table(paste(dirsel_boots_dir, "/boot_", i, ".bed", sep = ""), header = TRUE, sep = "\t")
    balsel_overl_file <- read.table(paste(balsel_boots_overl_dir, "/boot_", i, "_overl.bed", sep = ""), header = TRUE, sep = "\t")
    dirsel_overl_file <- read.table(paste(dirsel_boots_overl_dir, "/boot_", i, "_overl.bed", sep = ""), header = TRUE, sep = "\t")
    # Get total base pairs encompassed over the genome for bootstraps
    balsel_bp <- sum(unlist(apply(balsel_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
    dirsel_bp <- sum(unlist(apply(dirsel_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
    balsel_overl_bp <- sum(unlist(apply(balsel_overl_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))
    dirsel_overl_bp <- sum(unlist(apply(dirsel_overl_file, 1, function(row) as.numeric(row["end"]) - as.numeric(row["start"]))))

    # Get fraction of whole genome for each category
    fdiffs_frac <- fdiffs_bp/genome_size
    balsel_frac <- balsel_bp/genome_size
    dirsel_frac <- dirsel_bp/genome_size
    balsel_overl_frac <- balsel_overl_bp/genome_size
    dirsel_overl_frac <- dirsel_overl_bp/genome_size

    # add to storage vector
    balsel_overl_vec <- c(balsel_overl_vec, balsel_overl_bp)
    dirsel_overl_vec <- c(dirsel_overl_vec, dirsel_overl_bp)

    # Get expected balancing selection overlap and directional selection overlap expected values
    balsel_overl_expected <- fdiffs_frac * balsel_frac * genome_size
    dirsel_overl_expected <- fdiffs_frac * dirsel_frac * genome_size

    # add to storage vector
    balsel_overl_expected_vec <- c(balsel_overl_expected_vec, balsel_overl_expected)
    dirsel_overl_expected_vec <- c(dirsel_overl_expected_vec, dirsel_overl_expected)
}

# get difference between observed and expected overlap and the 95% confidence interval
balsel_diff_vec <- balsel_overl_vec - balsel_overl_expected_vec
dirsel_diff_vec <- dirsel_overl_vec - dirsel_overl_expected_vec

balsel_95conf_int <- sort(balsel_diff_vec)[c(3, 98)]
dirsel_95conf_int <- sort(dirsel_diff_vec)[c(3,98)]

# write results to disk
write.table(balsel_overl_vec, "boots/Eglob_all_balance_sel/Eglob_all_balance_sel_boots_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dirsel_overl_vec, "boots/Eglob_all_direction_sel/Eglob_all_direction_sel_boots_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(balsel_overl_expected_vec, "boots/Eglob_all_balance_sel/Eglob_all_balance_sel_boots_expected_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dirsel_overl_expected_vec, "boots/Eglob_all_direction_sel/Eglob_all_direction_sel_boots_expected_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(balsel_diff_vec, "Eglob_all_balance_boots_diff_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dirsel_diff_vec, "Eglob_all_direction_boots_diff_bp.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

95% of bootstrapped differences between observed and expected genome coverage fell within the following:

| Category           | boots diff lower bound (bp) | point estimate (bp) | boots diff upper bound (bp) |
| ------------------ | --------------------------- | ------------------- | --------------------------- |
| Balancing Sel      | -21333.378                  | -9715.5             | 7533.525                    |
| Directional Sel    | 106388.2                    | 167917.02           | 256741.5                    |

Which confirms that the balancing selection windows do not overlap more with species differences windows than would be expected by chance, while the directional selection windows do, by somewhere between ~100,000 and ~250,000 bp.

Plotted results in `R`.

```R
overlap_table <- data.frame(category = c("balancing", "directional"), 
                            lower_boot = c(-21333.378, 106388.2), 
                            point_est = c(-9715.5, 167917.02), 
                            upper_boot = c(7533.525, 256741.5))
# add margin space
par(mar = c(5, 4, 4, 8) + 0.1)
# plot point ests
overl_bar <- barplot(overlap_table$point_est, 
                     ylim = c(-50000, 300000),
                     col = c("#1ddda3", "#008a60"), 
                     main = "Overlap in selection and spp difference windows", 
                     ylab = "Overlap (bp)", 
                     cex.names = 1.2, 
                     legend.text = c("Balancing", "Directional", "95% bootstrap"),
                     args.legend = list(x = "topright", cex = 1.1, inset=c(-0.2, 0), 
                                        fill = c("#1ddda3", "#008a60", NA), 
                                        title = "Selection Type"))
# error bars
arrows(x0 = overl_bar, 
       y0 = overlap_table$upper_boot, 
       y1 = overlap_table$lower_boot, 
       angle = 90,
       code = 3,
       lwd = 2.5,
       length = 0.1)
# draw x axis line
abline(h = 0)
```

#### Average FST in selection windows
Bootstrapped combined selection window BED file in `R`.

```R
# import original BED files
fst_sel_MRcord_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/globMR_cordMR_fst_sel.bed"
fst_sel_purecord_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/synthesis/globpure_cordMR_fst_sel.bed"

fst_sel_MRcord_bed <- read.table(fst_sel_MRcord_filename, header = TRUE, sep = "\t")
fst_sel_purecord_bed <- read.table(fst_sel_purecord_filename, header = TRUE, sep = "\t")

# generate bootstrapped BED files
MRcord_boots_list <- list()
purecord_boots_list <- list()

for(i in c(1:100)){
    # generate and sort resampled indices
    fst_sel_MRcord_res_ind <- sort(sample(c(1:nrow(fst_sel_MRcord_bed)), nrow(fst_sel_MRcord_bed), replace = TRUE))
    fst_sel_purecord_res_ind <- sort(sample(c(1:nrow(fst_sel_purecord_bed)), nrow(fst_sel_purecord_bed), replace = TRUE))
    # subset original files with resampled indices
    fst_sel_MRcord_res <- fst_sel_MRcord_bed[fst_sel_MRcord_res_ind,] 
    fst_sel_purecord_res <- fst_sel_purecord_bed[fst_sel_purecord_res_ind,]
    # record to list
    MRcord_boots_list[[i]] <- fst_sel_MRcord_res
    purecord_boots_list[[i]] <- fst_sel_purecord_res
    # write resampled files
    write.table(fst_sel_MRcord_res, paste("boots/globMR_cordMR_fst_sel/boot_", i, ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    write.table(fst_sel_purecord_res, paste("boots/globpure_cordMR_fst_sel/boot_", i, ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

# Calculate average FST in bootstrapped BED files
# import FST table
fst_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/all_fst.txt"
fst_file <- read.table(fst_filename, header = TRUE, sep = "\t")

min_snps <- 15 # set minimum SNP number to consider a window

# define filters
# snp number filter
snp_fil <- which(fst_file$no_snps > min_snps)
# filter for correct population pairs
fst_MRcord_popfil <- which(fst_file$pop1 == "glob_MR" & fst_file$pop2 == "cord_MR")
fst_purecord_popfil <- which(fst_file$pop1 == "glob_pure" & fst_file$pop2 == "cord_MR")

# Loop through bootstraps, filtering FST file and calculating average FST for each
MRcord_fst_sel_avg_vec <- c()
purecord_fst_sel_avg_vec <- c()
for(i in c(1:100)){
    # get rows that are also overlapping with selection windows, as defined in imported BED file
    sel_MRcord_winfil <- unlist(apply(MRcord_boots_list[[i]], 1, function(row) which(fst_file$chromosome == row["chrom"] & fst_file$window_pos_1 == as.integer(row["start"]) & fst_file$window_pos_2 == as.integer(row["end"]))))
    sel_purecord_winfil <- unlist(apply(purecord_boots_list[[i]], 1, function(row) which(fst_file$chromosome == row["chrom"] & fst_file$window_pos_1 == as.integer(row["start"]) & fst_file$window_pos_2 == as.integer(row["end"]))))
    # final filters
    MRcord_final_fil <- intersect(snp_fil, intersect(fst_MRcord_popfil, sel_MRcord_winfil))
    purecord_final_fil <- intersect(snp_fil, intersect(fst_purecord_popfil, sel_purecord_winfil))

    # get final filtered tables for each population pair
    MRcord_fst_sel_tab <- fst_file[MRcord_final_fil,]
    purecord_fst_sel_tab <- fst_file[purecord_final_fil,]

    # store average FST
    MRcord_fst_sel_avg_vec <- c(MRcord_fst_sel_avg_vec, mean(MRcord_fst_sel_tab$avg_wc_fst))
    purecord_fst_sel_avg_vec <- c(purecord_fst_sel_avg_vec, mean(purecord_fst_sel_tab$avg_wc_fst))
}

# get final confidence intervals
MRcord_fst_sel_95conf_int <- sort(MRcord_fst_sel_avg_vec)[c(3,98)]
purecord_fst_sel_95conf_int <- sort(purecord_fst_sel_avg_vec)[c(3,98)]

# write bootstrapped averages to disk
write.table(MRcord_fst_sel_avg_vec, "MRcord_ft_sel_avg_boots.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(purecord_fst_sel_avg_vec, "purecord_ft_sel_avg_boots.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
```

Average FST values of bootstraps:
| Category            | boot lower bound | point estimate FST | boot upper bound | Genome-wide FST |
| ------------------- | ---------------- | ------------------ | ---------------- | --------------- |
| glob_MR & cord_MR   | 0.3504071        | 0.3658838          | 0.3805583        | 0.281874        |
| glob_pure & cord_MR | 0.3294071        | 0.3471521          | 0.3644682        | 0.2748824       |

95% of bootstrap average FSTs still fall well above the genome-wide FST value for both Meehan Range _E. globulus_ vs. _E. cordata_ and "pure" _E. globulus_ vs. _E. cordata_.

Plotted barplot of results in `R`.
```R
# make data frame of data
fst_table <- data.frame(category = c("globMR_cordMR", "globpure_cordMR"), 
                        lower_boot = c(0.3504071, 0.3294071), 
                        point_est = c(0.3658838, 0.3471521), 
                        upper_boot = c(0.3805583, 0.3644682), 
                        genome_wide = c(0.281874, 0.2748824))
# expand margins
par(mar = c(10, 4, 4, 2) + 0.1)
# plot barplot of point estimates
fst_bar <- barplot(fst_table$point_est, 
                   ylim = c(0, 0.5),
                   col = c("#7d7d7d", "#0087a8"), 
                   main = "Average FST in selection windows\n(balancing and directional combined)", 
                   ylab = "Average FST", 
                   cex.names = 1.2, 
                   legend.text = c("MR E. globulus - E. cordata", "pure E. globulus - E. cordata", "95% bootstrap", "Genome-wide Average"),
                   args.legend = list(x = "bottom", inset = c(0 , -0.325), 
                                      fill = c("#7d7d7d", "#0087a8", "NA", "NA"), 
                                      title = "Population Pair"))
# error bars
arrows(x0 = overl_bar, 
       y0 = fst_table$upper_boot, 
       y1 = fst_table$lower_boot, 
       angle = 90,
       code = 3,
       lwd = 2.5,
       length = 0.1)
# draw dashed lines for genome-wide averages
lines(x = c((fst_bar[1,1] - 0.5), (fst_bar[1,1] + 0.5)), 
      y = c(fst_table[1, "genome_wide"], fst_table[1, "genome_wide"]), 
      lty = 2, lwd = 2)
lines(x = c((fst_bar[2,1] - 0.5), (fst_bar[2,1] + 0.5)), 
      y = c(fst_table[2, "genome_wide"], fst_table[2, "genome_wide"]), 
      lty = 2, lwd = 2)
# draw line for x-axis
abline(h = 0)
```