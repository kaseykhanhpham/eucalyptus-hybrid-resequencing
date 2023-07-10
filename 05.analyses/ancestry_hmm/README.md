# Local Ancestry Inference
# HMM Inference of Local Ancestry
Performed using the tool [`Ancestry_HMM`](https://github.com/russcd/Ancestry_HMM). Thanks to Shelley Sianta for suggesting this analysis.

## Prepare inputs
### Genotype count input file

First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts.job for more details.
# Resources used: 5 Mb, 50 min

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out "$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### AncestryHMM input file
Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in.job for more details.
# Resources used: 4 Gb, 30 min

module load R/4.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt coverage.txt sample_files.txt gen_dists.tab all_ahmm_in.tab
```

## Run the program
Ran `Ancestry_HMM` using generated genotype counts input files and results from [`ADMIXTURE`](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/wg_ADMIXTURE) to provide contributions from source populations.

Initial run parameters:
* -i: input of genotype counts among reference and introgressed samples
* -s: ploidy of each sample, all diploid
* -a: two source populations, 0 (_E. globulus_) contributed 99% of variation and 1 (E. cordata) contributed 1%
* -p: first ancestry pulse for glob background, num generations set above limit to indicate starting background, 99% of current admixed genomes 
* -p: second ancestry pulse for cordy background, start at 1700 generations (assuming pleistocene admixture and generation times of 10 years) and optimize, 1% of current admixed genomes but optimize estimate
* -g: genotype counts provided rather than read pileups
* -b: do 100 bootstraps of 10,000 SNP blocks
* -ne: effective population size of the admixed population -- I don't know this, so I'm not going to try to provide it.

```bash
# Run in UFRC queue system; see ancestryhmm.job for more details.
# Resources used: 11 Gb, 3 hrs

module load ancestryhmm/1.0.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i "$WDIR"/all_ahmm_in.tab -s "$WDIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

## Analyze results

Plotted posteriors by chromosome for each admixed sample locally.
```bash
# Run in UFRC queue system; see plot_posteriors.job for more details.
# Resources used: 600 Mb, 10 min

module load R/4.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/plot_posteriors.r "$WDIR"/posteriors/"$NAME".posterior "goldenrod1,green3,deepskyblue4"
done < "$WDIR"/Eglobulus_MR.txt
```

### Overlap with other outlier windows

Retrieved windows with a posterior probability of > 95% for homozygous or heterozygous _E. cordata ancestry_.

```bash
# Run in UFRC queue system; see get_cord_posterior.job for more details.
# Resources used: 500 Mb, 7 min

module load R/4.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/get_cord_posterior.r "$WDIR"/posteriors/"$NAME".posterior "$NAME"_cord 0.95
done < "$WDIR"/Eglobulus_MR.txt
```

Converted posterior tables to genome annotation BED files.

```bash
module load R/4.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/post_to_bed.r "$NAME"_cord_het_0.95.tab bedfiles/"$NAME"_cord_het_0.95.bed_tmp
    Rscript "$SCRIPT_DIR"/post_to_bed.r "$NAME"_cord_hom_0.95.tab bedfiles/"$NAME"_cord_hom_0.95.bed_tmp
done < "$WDIR"/Eglobulus_MR.txt
```

Merged Ancestry_HMM loci together if <500 bp apart using `BEDTools`.

```bash
module load bedtools/2.30.0
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

while read NAME
do
    bedtools merge -d 500 -i "$NAME"_cord_het_0.95.bed_tmp > "$NAME"_cord_het_0.95.bed
    bedtools merge -d 500 -i "$NAME"_cord_hom_0.95.bed_tmp > "$NAME"_cord_hom_0.95.bed
done < "$LIST_DIR"/Eglobulus_MR.txt

rm *.bed_tmp
```

Examined overlap between `Ancestry_HMM`-inferred _cordata_ regions, regions of MR _globulus_ with low dxy, and regions with high fdM.

```bash
module load bedtools/2.30.0
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
DXY_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals/outlier_files/bedfiles"
AHMM_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/cord_loci/bedfiles"
DSUITE_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"

while read NAME
do
    # AHMM with DXY
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p95.bed > "$NAME"_ahmm_het_dxy_p95.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p90.bed > "$NAME"_ahmm_het_dxy_p90.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p95.bed > "$NAME"_ahmm_hom_dxy_p95.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p90.bed > "$NAME"_ahmm_hom_dxy_p90.bed

    # AHMM with f-stats
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DSUITE_FILES"/fDm_40_20_outliers_p05.bed > "$NAME"_ahmm_het_fDm.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DSUITE_FILES"/df_40_20_outliers_p05.bed > "$NAME"_ahmm_het_df.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DSUITE_FILES"/fDm_40_20_outliers_p05.bed > "$NAME"_ahmm_hom_fDm.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DSUITE_FILES"/df_40_20_outliers_p05.bed > "$NAME"_ahmm_hom_df.bed

done < "$LIST_DIR"/Eglobulus_MR.txt
```

### Characterize stats within high-posterior windows

Converted `pixy` output into BED file of windows.

```R
wdir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
setwd(wdir)
options(scipen=999)

all_dxy <- read.table("all_dxy.txt", header = TRUE)
all_pi <- read.table("all_pi.txt", header = TRUE)

dxy_bed <- all_dxy[which(all_dxy$pop1 == "glob_MR" & all_dxy$pop2 == "cord_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
dxy_bed$window_pos_1 <- as.numeric(dxy_bed$window_pos_1) - 1
write.table(dxy_bed, "all_dxy.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

pi_bed <- all_pi[which(all_pi$pop == "glob_MR"), c("chromosome", "window_pos_1", "window_pos_2")]
pi_bed$window_pos_1 <- as.numeric(pi_bed$window_pos_1) - 1
write.table(pi_bed, "all_pi.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
```

Got pi/dxy windows for windows of high posterior probability of _cordata_ ancestry.

```bash
module load bedtools/2.30.0
STAT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
AHMM_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/cord_loci/bedfiles"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

while read NAME
do
    bedtools intersect -a "$STAT_DIR"/all_pi.bed -b "$AHMM_DIR"/"$NAME"_cord_hom_0.95.bed -wa > "$NAME"_cord_hom_0.95_pixy_windows.bed
    bedtools intersect -a "$STAT_DIR"/all_pi.bed -b "$AHMM_DIR"/"$NAME"_cord_het_0.95.bed -wa > "$NAME"_cord_het_0.95_pixy_windows.bed
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Retrieved pi/dxy values for windows.

```bash
module load R/4.2
PIXY_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
BED_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/post_stats/pixy_windows"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

# summary stats for E. glob MR and all E. cordata (looking for high values of pi and low values of dxy/fst)
# Minimum of 40 SNPs per window
touch ahmm_window_stats_40.txt
while read NAME
do
    # pi
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_pi.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 40 "pi" "glob_MR" >> ahmm_window_stats_40.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_pi.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 40 "pi" "glob_MR">> ahmm_window_stats_40.txt

    # dxy
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 40 "dxy" "glob_MR" "cord_MR" >> ahmm_window_stats_40.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 40 "dxy" "glob_MR" "cord_MR" >> ahmm_window_stats_40.txt

    # fst
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_fst.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 40 "fst" "glob_MR" "cord_MR" >> ahmm_window_stats_40.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_fst.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 40 "fst" "glob_MR" "cord_MR" >> ahmm_window_stats_40.txt
done < "$LIST_DIR"/Eglobulus_MR.txt

# Minimum of 20 SNPs per winndow (little change from 40 SNPs per window)
touch ahmm_window_stats_20.txt
while read NAME
do
    # pi
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_pi.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 20 "pi" "glob_MR" >> ahmm_window_stats_20.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_pi.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 20 "pi" "glob_MR">> ahmm_window_stats_20.txt

    # dxy
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 20 "dxy" "glob_MR" "cord_MR" >> ahmm_window_stats_20.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 20 "dxy" "glob_MR" "cord_MR" >> ahmm_window_stats_20.txt

    # fst
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_fst.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 20 "fst" "glob_MR" "cord_MR" >> ahmm_window_stats_20.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_fst.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 20 "fst" "glob_MR" "cord_MR" >> ahmm_window_stats_20.txt
done < "$LIST_DIR"/Eglobulus_MR.txt

# dxy between individuals in E. glob MR and all E. cordata (looking for low values)
PIXY_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals"
while read NAME
do
    # dxy
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/cord_"$NAME"_dxy.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 40 "dxy" "glob_MR" "cord_MR" >> ahmm_window_stats_inds.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/cord_"$NAME"_dxy.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 40 "dxy" "cord_MR" "glob_MR" >> ahmm_window_stats_inds.txt
done < "$LIST_DIR"/Eglobulus_MR.txt

# dxy between E. glob MR and E. glob reference (looking for high values)
PIXY_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
while read NAME
do
    # dxy
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_hom_0.95_pixy_windows.bed 40 "dxy" "glob_MR" "glob_pure" >> ahmm_window_stats_glob.txt
    Rscript "$SCRIPT_DIR"/calc_window_stats.r "$PIXY_DIR"/all_dxy.txt "$BED_DIR"/"$NAME"_cord_het_0.95_pixy_windows.bed 40 "dxy" "glob_MR" "glob_pure" >> ahmm_window_stats_glob.txt
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Most windows with pi/dxy at or a bit larger than genome-wide calculation... Could be that introgressed alleles are at a low enough frequency in _E. glob_ MR population to not show signal in dxy. 

### Visualize shared loci

Split chromosomes into 5 sections each and plotted presence/absence of heterozygous and homozygous _cordata_ loci above 95% posterior probability.

```R
library(ggplot2)
library(RColorBrewer)

post_file_loc <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/posteriors"
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/heatmaps"
glob_mr_samples <- c("WA01", "WA03", "WA04", "WB02", "WB03", "WB04", "WC02", "WC03", "WC05", "WD04", "WE02", "WE03", "WE04", "WE05", "WF01", "WG03", "WG04", "WG05", "WH03", "WH04")
chromosomes <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

populate_mat <- function(sample_vec, window_tab){
    mat_list <- list()
    index = 1
    for(sample in sample_vec){
        infile_name <- paste(post_file_loc, "/", sample, ".posterior", sep = "")
        infile_colnames <- c("chrom", "pos", "hom_glob", "het", "hom_cord")
        infile <- read.table(infile_name, skip = 1, sep = "\t", col.names = infile_colnames, colClasses = c("character", "integer", "numeric", "numeric", "numeric"))
        # construct row for each window -- 1 for a cell if the ind has any homozygous cord posterior > 95% in the window
        # 0.5 for a cell if the ind has any heterozygous cord posterior > 95% in the window
        for(i in c(1:nrow(window_tab))) {
            wchr <- window_tab[i, "chrom"]
            wsta <- window_tab[i, "start"]
            wend <- window_tab[i, "end"]
            infile_mask <- which(infile$chrom == wchr & infile$pos >= wsta & infile$pos <= wend)
            infile_subset <- infile[infile_mask,]
            if(any(infile_subset$hom_cord >= 0.95)){
                hom_cell <- 1
                het_cell <- 1
            } else if(any(infile_subset$het >= 0.95)){
                hom_cell <- 0
                het_cell <- 0.5
            } else {
                hom_cell <- 0
                het_cell <- 0
            }
            windex <- paste(wchr, wsta, wend, sep = "_")
            mat_ind <- paste("row", index, sep = "_")

            mat_list[[mat_ind]]["window"] <- windex
            mat_list[[mat_ind]]["sample"] <- sample
            mat_list[[mat_ind]]["hom_mat"] <- hom_cell
            mat_list[[mat_ind]]["het_mat"] <- het_cell
            
            index <- index + 1
        }
    }
    mat_df <- as.data.frame(t(as.data.frame(mat_list)))
    mat_df$hom_mat <- as.integer(mat_df$hom_mat)
    mat_df$het_mat <- as.numeric(mat_df$het_mat)
    return(mat_df)
}

setwd(wdir)
coarsewin_filename_rel <- "coarse_windows.csv"
coarsewin_filename <- paste(wdir, window_filename_rel, sep = "/")
coarsewin_tab <- read.csv(window_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))

# populate matrix
coarse_mat_df <- populate_mat(glob_mr_samples, coarsewin_tab)

# plot heatmaps
x_axis_breaks <- c("Chr01_16887823_25331733", "Chr02_20331353_30497028", "Chr03_26218899_39328347", "Chr04_15439735_23159601"
, "Chr05_25167991_37751985", "Chr06_20856265_31284396", "Chr07_21701053_32551578", "Chr08_28085845_42128766", "Chr09_15320131_22980195", "Chr10_15489065_23233596", "Chr11_16822585_25233876")

name_table_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
name_table <- read.csv(name_table_name, header = TRUE, as.is = TRUE)
label_order <- match(coarse_mat_df$sample, name_table$RAPiD_ID)
coarse_mat_df$acc <- name_table$Accession[label_order]

hom_coarse <- ggplot(mat_df, aes(window, acc, fill = hom_mat)) + 
              geom_tile(color = "#00617a") + 
              scale_fill_steps(high = "#ffd966", low = "#0a9cc1") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
hom_coarse

het_coarse <- ggplot(mat_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
het_coarse
```

![Presence of homozygous _E. cordata_ ancestry per chromosome and sample; chromosome 6 is a hotspot.](coarse_hom_windows.png "Homozygous _E. cordata_ ancestry heatmap")

![Presence of heterozygous or homozygous _E. cordata_ ancestry per chromosome and sample; heterozygotes were scored as half in presence/absence matrix.](coarse_het_windows.png "_E. cordata_ ancestry heatmap")

Blocks of interest for homozygous _E. cordata_ ancestry at the following based on sharing between individuals:

| Chromosome | Start (bp)  | End (bp)    |
| ---------- | ----------- | ----------- |
| Chr06      | 10,428,133  | 20,856,264  |
| Chr07      | 1           | 10,850,526  |
| Chr08      | 14,042,923  | 28,085,844  |

Split blocks of interest into finer windows.

```R
library(ggplot2)

hm_theme <- theme(plot.title=element_text(size=16),
                  axis.title = element_text(size = 14),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 14))

# Chromosome 6, Block 2
chr06b2_samples <- c("WG05", "WC03", "WC02", "WB04", "WE04", "WH03")
chr06b2_filename_rel <- "chr06_block2_windows.csv"
chr06b2_filename <- paste(wdir, chr06b2_filename_rel, sep = "/")
chr06b2_tab <- read.csv(chr06b2_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr06b2_df <- populate_mat(chr06b2_samples, chr06b2_tab)

label_order <- match(chr06b2_df$sample, name_table$RAPiD_ID)
chr06b2_df$acc <- name_table$Accession[label_order]
chr06b2_xlabs <- unique(gsub("Chr06_", "", chr06b2_df$window))

chr06b2_hom <- ggplot(chr06b2_df, aes(window, acc, fill = hom_mat)) + 
              geom_tile(color = "#00617a") + 
              scale_fill_steps(high = "#ffd966", low = "#0a9cc1") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry") + 
              xlab("Chromosome 06 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr06b2_xlabs) +
              hm_theme
chr06b2_hom
chr06b2_het <- ggplot(chr06b2_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 06 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr06b2_xlabs) +
              hm_theme
chr06b2_het

# Chromosome 7, Block 1
chr07b1_samples <- c("WB02", "WH04", "WG03", "WC03", "WE02", "WC02", "WB04", "WE04", "WA01", "WA04", "WH03", "WG04")
chr07b1_filename_rel <- "chr07_block1_windows.csv"
chr07b1_filename <- paste(wdir, chr07b1_filename_rel, sep = "/")
chr07b1_tab <- read.csv(chr07b1_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_df <- populate_mat(chr07b1_samples, chr07b1_tab)

label_order <- match(chr07b1_df$sample, name_table$RAPiD_ID)
chr07b1_df$acc <- name_table$Accession[label_order]
chr07b1_xlabs <- unique(gsub("Chr07_", "", chr07b1_df$window))

chr07b1_hom <- ggplot(chr07b1_df, aes(window, acc, fill = hom_mat)) + 
              geom_tile(color = "#00617a") + 
              scale_fill_steps(high = "#ffd966", low = "#0a9cc1") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_xlabs) +
              hm_theme
chr07b1_hom
chr07b1_het <- ggplot(chr07b1_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_xlabs) +
              hm_theme
chr07b1_het

# Chromosome 8, Block 2
chr08b2_samples <- c("WB02", "WG05", "WE02", "WC02", "WE04")
chr08b2_filename_rel <- "chr08_block2_windows.csv"
chr08b2_filename <- paste(wdir, chr08b2_filename_rel, sep = "/")
chr08b2_tab <- read.csv(chr08b2_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b2_df <- populate_mat(chr08b2_samples, chr08b2_tab)

label_order <- match(chr08b2_df$sample, name_table$RAPiD_ID)
chr08b2_df$acc <- name_table$Accession[label_order]
chr08b2_xlabs <- unique(gsub("Chr08_", "", chr08b2_df$window))

chr08b2_hom <- ggplot(chr08b2_df, aes(window, acc, fill = hom_mat)) + 
              geom_tile(color = "#00617a") + 
              scale_fill_steps(high = "#ffd966", low = "#0a9cc1") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry") + 
              xlab("Chromosome 08 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr08b2_xlabs) +
              hm_theme
chr08b2_hom
chr08b2_het <- ggplot(chr08b2_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 08 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr08b2_xlabs) +
              hm_theme
chr08b2_het
```

![Presence of _E. cordata_ ancestry in Chromosome 6 10428133 bp to 20856264 bp; heterozygotes were scored as half in presence/absence matrix.](chr06b2_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 6 Block 2")

![Presence of _E. cordata_ ancestry in Chromosome 7 1 bp to 10850526 bp; heterozygotes were scored as half in presence/absence matrix.](chr07b1_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1")

![Presence of _E. cordata_ ancestry in Chromosome 8 14042923 to 28085844 bp; heterozygotes were scored as half in presence/absence matrix.](chr08b2_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 8 Block 2")

Blocks of interest for homozygous _E. cordata_ ancestry at the following based on sharing between individuals:

| Chromosome | Start (bp)  | End (bp)    |
| ---------- | ----------- | ----------- |
| Chr06      | 12,513,761  | 13,556,574  |
| Chr07      | 1           | 1,085,053   |
| Chr07      | 6,510,319   | 7,595,371   |
| Chr08      | 21,064,388  | 23,872,973  |

Plotted zoomed shared blocks, split again.

```R
library(ggplot2)

# Chromosome 6 Block 2-3
chr06b2_3_filename_rel <- "chr06_block2-3_windows.csv"
chr06b2_3_filename <- paste(wdir, chr06b2_3_filename_rel, sep = "/")
chr06b2_3_tab <- read.csv(chr06b2_3_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr06b2_3_df <- populate_mat(chr06b2_samples, chr06b2_3_tab)

label_order <- match(chr06b2_3_df$sample, name_table$RAPiD_ID)
chr06b2_3_df$acc <- name_table$Accession[label_order]
chr06b2_3_xlabs <- unique(gsub("Chr06_", "", chr06b2_3_df$window))

chr06b2_3_het <- ggplot(chr06b2_3_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 06 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr06b2_3_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr06b2_3_het

# Chromosome 7 Block 1-1
chr07b1_1_filename_rel <- "chr07_block1-1_windows.csv"
chr07b1_1_filename <- paste(wdir, chr07b1_1_filename_rel, sep = "/")
chr07b1_1_tab <- read.csv(chr07b1_1_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_1_df <- populate_mat(chr07b1_samples, chr07b1_1_tab)

label_order <- match(chr07b1_1_df$sample, name_table$RAPiD_ID)
chr07b1_1_df$acc <- name_table$Accession[label_order]
chr07b1_1_xlabs <- unique(gsub("Chr07_", "", chr07b1_1_df$window))

chr07b1_1_het <- ggplot(chr07b1_1_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_1_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr07b1_1_het

# Chromosome 7 Block 1-7
chr07b1_7_filename_rel <- "chr07_block1-7_windows.csv"
chr07b1_7_filename <- paste(wdir, chr07b1_7_filename_rel, sep = "/")
chr07b1_7_tab <- read.csv(chr07b1_7_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_7_df <- populate_mat(chr07b1_samples, chr07b1_7_tab)

label_order <- match(chr07b1_7_df$sample, name_table$RAPiD_ID)
chr07b1_7_df$acc <- name_table$Accession[label_order]
chr07b1_7_xlabs <- unique(gsub("Chr07_", "", chr07b1_7_df$window))

chr07b1_7_het <- ggplot(chr07b1_7_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_7_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr07b1_7_het

# Chromosome 8 Block 2-67
chr08b2_67_filename_rel <- "chr08_block2-67_windows.csv"
chr08b2_67_filename <- paste(wdir, chr08b2_67_filename_rel, sep = "/")
chr08b2_67_tab <- read.csv(chr08b2_67_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b2_67_df <- populate_mat(chr08b2_samples, chr08b2_67_tab)

label_order <- match(chr08b2_67_df$sample, name_table$RAPiD_ID)
chr08b2_67_df$acc <- name_table$Accession[label_order]
chr08b2_67_xlabs <- unique(gsub("Chr08_", "", chr08b2_67_df$window))

chr08b2_67_het <- ggplot(chr08b2_67_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 08 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr08b2_67_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr08b2_67_het
```

![Presence of _E. cordata_ ancestry in Chromosome 6 12513761 bp to 13556574 bp; heterozygotes were scored as half in presence/absence matrix.](chr06b2-3_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 6 Block 2_3")

![Presence of _E. cordata_ ancestry in Chromosome 7 1 bp to 1085053 bp; heterozygotes were scored as half in presence/absence matrix.](chr07b1-1_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1_1")

![Presence of _E. cordata_ ancestry in Chromosome 7 6510319 bp to 7595371 bp; heterozygotes were scored as half in presence/absence matrix.](chr07b1-7_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1_7")

![Presence of _E. cordata_ ancestry in Chromosome 8 21064388 to 23872973 bp; heterozygotes were scored as half in presence/absence matrix.](chr08b2-67_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 8 Block 2_6-7")

| Chromosome | Start (bp)  | End (bp)    |
| ---------- | ----------- | ----------- |
| Chr06      | 13,139,453  | 13,190,593  |
| Chr07      | 217,013     | 488,277     |
| Chr07      | 7,052,849   | 7,161,354   |
| Chr08      | 22,328,249  | 22,749,536  |

Plotted zoomed shared blocks, split again.

```R
library(ggplot2)

# Chromosome 6 Block 2-3
chr06b2_3_13_filename_rel <- "chr06_block2-3-13_windows.csv"
chr06b2_3_13_filename <- paste(wdir, chr06b2_3_13_filename_rel, sep = "/")
chr06b2_3_13_tab <- read.csv(chr06b2_3_13_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr06b2_3_13_df <- populate_mat(chr06b2_samples, chr06b2_3_13_tab)

label_order <- match(chr06b2_3_13_df$sample, name_table$RAPiD_ID)
chr06b2_3_13_df$acc <- name_table$Accession[label_order]
chr06b2_3_13_xlabs <- unique(gsub("Chr06_", "", chr06b2_3_13_df$window))

chr06b2_3_13_het <- ggplot(chr06b2_3_13_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 06 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr06b2_3_13_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr06b2_3_13_het

# Chromosome 7 Block 1-1
chr07b1_1_59_filename_rel <- "chr07_block1-1-59_windows.csv"
chr07b1_1_59_filename <- paste(wdir, chr07b1_1_59_filename_rel, sep = "/")
chr07b1_1_59_tab <- read.csv(chr07b1_1_59_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_1_59_df <- populate_mat(chr07b1_samples, chr07b1_1_59_tab)

label_order <- match(chr07b1_1_59_df$sample, name_table$RAPiD_ID)
chr07b1_1_59_df$acc <- name_table$Accession[label_order]
chr07b1_1_59_xlabs <- unique(gsub("Chr07_", "", chr07b1_1_59_df$window))

chr07b1_1_59_het <- ggplot(chr07b1_1_59_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_1_59_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr07b1_1_59_het

# Chromosome 7 Block 1-7
chr07b1_7_1112_filename_rel <- "chr07_block1-7-1112_windows.csv"
chr07b1_7_1112_filename <- paste(wdir, chr07b1_7_1112_filename_rel, sep = "/")
chr07b1_7_1112_tab <- read.csv(chr07b1_7_1112_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_7_1112_df <- populate_mat(chr07b1_samples, chr07b1_7_1112_tab)

label_order <- match(chr07b1_7_1112_df$sample, name_table$RAPiD_ID)
chr07b1_7_1112_df$acc <- name_table$Accession[label_order]
chr07b1_7_1112_xlabs <- unique(gsub("Chr07_", "", chr07b1_7_1112_df$window))

chr07b1_7_1112_het <- ggplot(chr07b1_7_1112_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 07 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr07b1_7_1112_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr07b1_7_1112_het

# Chromosome 8 Block 2-67
chr08b2_67_1012_filename_rel <- "chr08_block2-67-1012_windows.csv"
chr08b2_67_1012_filename <- paste(wdir, chr08b2_67_1012_filename_rel, sep = "/")
chr08b2_67_1012_tab <- read.csv(chr08b2_67_1012_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b2_67_1012_df <- populate_mat(chr08b2_samples, chr08b2_67_1012_tab)

label_order <- match(chr08b2_67_1012_df$sample, name_table$RAPiD_ID)
chr08b2_67_1012_df$acc <- name_table$Accession[label_order]
chr08b2_67_1012_xlabs <- unique(gsub("Chr08_", "", chr08b2_67_1012_df$window))

chr08b2_67_1012_het <- ggplot(chr08b2_67_1012_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry") + 
              xlab("Chromosome 08 Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(labels = chr08b2_67_1012_xlabs) +
              hm_theme + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
chr08b2_67_1012_het
```

![Presence of _E. cordata_ ancestry in Chromosome 6 13139453 bp to 13190593 bp; heterozygotes were scored as half in presence/absence matrix.](chr06b2-3-13_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 6 Block 2_3_13")

![Presence of _E. cordata_ ancestry in Chromosome 7 217013 bp to 488277 bp; heterozygotes were scored as half in presence/absence matrix.](chr07b1-1-59_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1_1_5-9")

![Presence of _E. cordata_ ancestry in Chromosome 7 7052849 bp to 7161354 bp; heterozygotes were scored as half in presence/absence matrix.](chr07b1-7-1112_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1_7_11-12")

![Presence of _E. cordata_ ancestry in Chromosome 8 22328249 to 22749536 bp; heterozygotes were scored as half in presence/absence matrix.](chr08b2-67-1012_het_windows.png "_E. cordata_ ancestry heatmap for Chromosome 8 Block 2_6-7_10-12")

Final intervals to investigate further:

| Chromosome | Start (bp)  | End (bp)    |
| ---------- | ----------- | ----------- |
| Chr06      | 13,159,909  | 13,180,364  |
| Chr07      | 230,577     | 237,358     |
| Chr07      | 278,051     | 284,832     |
| Chr07      | 427,255     | 447,600     |
| Chr07      | 7,074,553   | 7,123,386   |
| Chr08      | 22,454,645  | 22,486,243  |

Intervals for single samples:

| Sample | Chromosome | Start (bp)  | End (bp)    |
| ------ | ---------- | ----------- | ----------- |
| 6024   | Chr08      | 22,496,777  | 22,602,106  |
| 4207   | Chr08      | 22,654,772  | 22,739,035  |
