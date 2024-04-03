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
# Set population pair to to correct populations

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

| Category                          | bp covered      | fraction of genome | expected bp covered |
| --------------------------------- | --------------- | ------------------ | ------------------- |
| species differences               | 5958976         | 0.009877278        | NA                  |
| balancing selection               | 2501955         | 0.004147106        | NA                  |
| directional selection             | 3042942         | 0.005043817        | NA                  |
| spp diffs & balancing selection   | 14997           | 0.00002485822      | 24712.5             |
| spp diffs & directional selection | 197973          | 0.0003281494       | 30055.98            |

The amount of balancing selection windows overlapping with species difference windows is about what you'd expect from random chance. However, the amount of directional selection windows overlapping with species difference windows is much higher, by a factor of ~6. This indicates that divergent selection on _E. globulus_ may be acting against any introgression from _E. cordata_.

As for average FST of selection windows...
| Category            | FST           | Genome-wide FST |
| ------------------- | ------------- | --------------- |
| glob_MR & cord_MR   | 0.3658838     | 0.281874        |
| glob_pure & cord_MR | 0.3471521     | 0.2748824       |

Average FST in both types of selection windows combined is also quite a bit higher than genome-wide FST, filtering for a minimum number of 15 SNPs per window (which of course biases the numbers towards regions of the genome with more coverage).

### Error Estimation
Bootstrapped selection windows to get errors on point estimates of species difference/selection window overlap and average FST in selection windows. Started by generating the bootstraps in `R`.

```R
## Compare overlap between species difference and selection windows ##
#  ----------------------------------------------------------------- #
# Import BED files as data frames
balsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/balancing/Eglob_all_balance_sel.bed"
dirsel_filename <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/sel_files/directional/Eglob_all_direction_sel.bed"

balsel_file <- read.table(balsel_filename, header = TRUE, sep = "\t")
dirsel_file <- read.table(dirsel_filename, header = TRUE, sep = "\t")

# Resample selection windows (100x)
balsel_boots <- list()
dirsel_boots <- list()

for(i in c(1:100)){
    balsel_boots[[i]] <- sample(c(1:nrow(balsel_file)), 100, replace=TRUE)
    dirsel_boots[[i]] <- sample(c(1:nrow(dirsel_file)), 100, replace = TRUE)
}
```