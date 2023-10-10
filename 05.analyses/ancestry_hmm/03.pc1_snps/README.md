# HMM Local Ancestry Inference on just PC-Correlated SNPs

## Data Preparation
### Filtering and Genotype Counts
Exported positions of PC1-correlated SNPs from `PCAdapt` analysis present in all samples as a tab-separated table of chromosome and physical position. Used this table to subset during genotype counts.

```bash
# Run in UFRC queue system; see gt_counts_pc.job for more details.
# Resources used: 

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --positions "$POPLIST_DIR"/03.pc1_snps/pc1_snps_pos.tab --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --positions "$POPLIST_DIR"/03.pc1_snps/pc1_snps_pos.tab --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out "$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --positions "$POPLIST_DIR"/03.pc1_snps/pc1_snps_pos.tab --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### Format AHMM Input Files
Used custom python script to format inputs for AHMM.

```bash
# Run in UFRC queue system; see get_ahmm_in_pc.job for more details.
# Resources used:

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/03.pc1_snps/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt "$LIST_DIR"/03.pc1_snps/coverage.txt sample_files.txt gen_dists.tab all_ahmm_in_pc.tab
```

## Run Ancestry_HMM
```bash
# Run in UFRC queue system; see ancestryhmm_pc.job for more details.
# Resources used:

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
ancestry_hmm -i all_ahmm_in_pc.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

## Summary Plots
### Across all chromosomes
Split chromosomes into five chunks each and marked any windows with at least one variant with posterior > 0.95 of homozygous _E. cordata_ ancestry.

```R
library(ggplot2)
library(RColorBrewer)

source("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/heatmaps_fun.r")

post_file_loc <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/03.pc1_snps/posteriors"
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/03.pc1_snps/heatmaps"
glob_mr_samples <- c("WA01", "WA03", "WA04", "WB02", "WB03", "WB04", "WC02", "WC03", "WC05", "WD04", "WE02", "WE03", "WE04", "WE05", "WF01", "WG03", "WG04", "WG05", "WH03", "WH04")
chromosomes <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

setwd(wdir)
coarsewin_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/high_cov/heatmaps/coarse_windows.csv"
coarsewin_tab <- read.csv(coarsewin_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))

# populate matrix
coarse_mat_df <- populate_mat(glob_mr_samples, coarsewin_tab, 0.95, post_file_loc)

# plot heatmaps
x_axis_breaks <- c("Chr01_16887823_25331733", "Chr02_20331353_30497028", "Chr03_26218899_39328347", "Chr04_15439735_23159601"
, "Chr05_25167991_37751985", "Chr06_20856265_31284396", "Chr07_21701053_32551578", "Chr08_28085845_42128766", "Chr09_15320131_22980195", "Chr10_15489065_23233596", "Chr11_16822585_25233876")

name_table_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
name_table <- read.csv(name_table_name, header = TRUE, as.is = TRUE)
label_order <- match(coarse_mat_df$sample, name_table$RAPiD_ID)
coarse_mat_df$acc <- name_table$Accession[label_order]

hom_coarse <- ggplot(coarse_mat_df, aes(window, acc, fill = hom_mat)) + 
              geom_tile(color = "#00617a") + 
              scale_fill_steps(high = "#ffd966", low = "#0a9cc1") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry - PC1 SNPs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
hom_coarse

het_coarse <- ggplot(coarse_mat_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry - PC1 SNPs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
het_coarse
```

![Presence of homozygous _E. cordata_ ancestry > 0.95 posterior per chromosome and sample; chromosomes 1,8 have regions shared between a large number of individuals.](heatmaps/coarse_hom_windows.png "Homozygous _E. cordata_ ancestry heatmap, PC1-correlated SNPs")

Results look completely different from using all variants or the AIMs. The following were regions of interest:

| Chromosome | Start (bp)  | End (bp)    | Number of Inds Sharing |
| ---------- | ----------- | ----------- | ---------------------- |
| Chr01      | 33,775,645  | 42,219,553  | 7                      |
| Chr08      | 28,085,845  | 42,128,766  | 6                      |

### Finer Windows
```R
# Chromosome 1, Block 5
chr01b5_samples <- c("WC05", "WE02", "WC02", "WB03", "WH03", "WG04", "WE05")
chr01b5_filename_rel <- "chr01_b5.csv"
chr01b5_filename <- paste(wdir, chr01b5_filename_rel, sep = "/")
chr01b5_tab <- read.csv(chr01b5_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr01b5_df <- populate_mat(chr01b5_samples, chr01b5_tab, 0.95, post_file_loc)

chr01b5_plots <- plot_heatmap(chr01b5_df, name_table$RAPiD_ID, name_table$Accession, "Chr01", 0.95, "- Chr01 B5 - PC1 SNPs")
chr01b5_plots$het_plot

# Chromosome 8, Block 3
chr08b3_samples <- c("WH04", "WG03", "WC03", "WC02", "WE04", "WE05")
chr08b3_filename_rel <- "chr08_b3.csv"
chr08b3_filename <- paste(wdir, chr08b3_filename_rel, sep = "/")
chr08b3_tab <- read.csv(chr08b3_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b3_df <- populate_mat(chr08b3_samples, chr08b3_tab, 0.95, post_file_loc)

chr08b3_plots <- plot_heatmap(chr08b3_df, name_table$RAPiD_ID, name_table$Accession, "Chr08", 0.95, "- Chr08 B3 - AIMs")
chr08b3_plots$het_plot
```

Split into even finer blocks:

```R
# Chromosome 1 Block 5_7
chr01b5_7_filename_rel <- "chr01_b5_7.csv"
chr01b5_7_filename <- paste(wdir, chr01b5_7_filename_rel, sep = "/")
chr01b5_7_tab <- read.csv(chr01b5_7_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr01b5_7_df <- populate_mat(chr01b5_samples, chr01b5_7_tab, 0.95, post_file_loc)

chr01b5_7_plots <- plot_heatmap(chr01b5_7_df, name_table$RAPiD_ID, name_table$Accession, "Chr01", 0.95, "- Chr01 B5 - PC1 SNPs")
chr01b5_7_plots$het_plot

# Chromosome 8 Block 3_1
chr08b3_1_filename_rel <- "chr08_b3_1.csv"
chr08b3_1_filename <- paste(wdir, chr08b3_1_filename_rel, sep = "/")
chr08b3_1_tab <- read.csv(chr08b3_1_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b3_1_df <- populate_mat(chr08b3_samples, chr08b3_1_tab, 0.95, post_file_loc)

chr08b3_1_plots <- plot_heatmap(chr08b3_1_df, name_table$RAPiD_ID, name_table$Accession, "Chr08", 0.95, "- Chr08 B3 - AIMs")
chr08b3_1_plots$het_plot
```

Once zoomed in, only three individuals are homozygous at a single block for the _E. cordata_ variants for both windows of interest, so I am stopping for now.