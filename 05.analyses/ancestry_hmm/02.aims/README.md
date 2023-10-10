# HMM Inference of Local Ancestry - Ancestry-Informative Markers
Several studies using `Ancestry_HMM` filtered SNPs first by informativeness in ancestral populations, either using FST or allele frequency. Since `AHMM` uses allele frequency to assign each individual's loci to a population, this should increase the power of the study and my confidence that we're not inferring off of high-conservation regions.

## High-Coverage Analysis
### Data Preparation
#### Marker Filtering

Plotted average Weir-Cockerham FST values [calculated in `pixy`](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/genome_scan/pixy) from 5kb windows across the genome between the _E. globulus_ reference group and _E. cordata_ to determine filtering cutoff. FST values appear in a right-skewed distribution with a median of 0.2641 and a mean of 0.2934. As these species are long-diverged and admixture is expected to have occurred only intermittently, I decided that a cutoff at the genome-wide median would be fine for capturing SNPs which would differentiate between species.

![Histogram of average Weir-Cockerham FST values from 5kb windows across genome between reference E. globulus group and E. cordata. The distribution goes from -0.1232 to 1.00 with a light right skew and a mean at 0.2934. The histogram meets its maximum height about 6000 windows just before the median and mean.](fst_sourcepops_distr.png "Average window FST between reference E. globulus and E. cordata")

Extracted windows above the median FST in `R` to BED format.

```R
fst_tab_name <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/all_fst.txt"
fst_tab <- read.table(fst_tab_name, header = TRUE)

aims_filter <- which(fst_tab$pop1 == "glob_pure" & fst_tab$pop2 == "cord_MR" & fst_tab$no_snps > 15 & fst_tab$avg_wc_fst > 0.2641)
fst_bed_fil <- fst_tab[aims_filter, c("chromosome", "window_pos_1", "window_pos_2")]
fst_bed_fil$window_pos_1 <- fst_bed_fil$window_pos_1 - 1

fst_bed_merged <- data.frame(chromosome = character(), window_pos_1 = integer(), window_pos_2 = integer())
for(i in c(nrow(fst_bed_fil):2)){
    if(fst_bed_fil[i,"window_pos_1"] == fst_bed_fil[(i-1),"window_pos_2"]){
        fst_bed_fil[(i-1),"window_pos_2"] <- fst_bed_fil[i, "window_pos_2"]
    } else {
        fst_bed_merged[(nrow(fst_bed_merged) + 1),] <- fst_bed_fil[i,]
    }
}
fst_bed_merged[(nrow(fst_bed_merged) + 1),] <- fst_bed_fil[1,]
fst_bed_merged_ordered <- fst_bed_merged[c(nrow(fst_bed_merged):1),]

colnames(fst_bed_merged_ordered) <- c("chrom", "chromStart", "chromEnd")
write.table(fst_bed_merged_ordered, "filtered_aims.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```

### AHMM Data Formatting
Got allele counts for all individuals in source populations and focal population, filtering sites by criteria above.

```bash
# Run in UFRC queue system; see gt_counts_aims.job for more details.
# Resources used: 

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --bed filtered_aims.bed --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --bed filtered_aims.bed --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out "$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --bed filtered_aims.bed --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

Generated input for running Ancestry_HMM using a custom python script.

```bash
# Run in UFRC queue system; see get_ahmm_in_aims.job for more details.
# Resources used: 

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/02.aims/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt "$LIST_DIR"/coverage.txt sample_files.txt gen_dists.tab all_ahmm_in_aims.tab
```

### Run Program and Visualize Results
#### Run Ancestry_HMM

Arguments:
* -i: input of genotype counts among reference and introgressed samples
* -s: ploidy of each sample, all diploid
* -a: two source populations, 0 (E. globulus) contributed 99% of variation and 1 (E. cordata) contributed 1%
* -p: first ancestry pulse for glob background, num generations set above limit to indicate starting background, 99% of current admixed genomes 
* -p: second ancestry pulse for cordy background, start at 1700 generations (assuming pleistocene admixture and generation times of 10 years) and optimize, 1% of current admixed genomes but optimize estimate
* -g: genotype counts provided rather than read pileups
* -b: do 100 bootstraps of 10,000 SNP blocks
* -ne: effective population size of the admixed population -- I don't know this, so I'm not going to try to provide it.

```bash
# Run in UFRC queue system; see ancestry_hmm_aims.job for more details.
# Resources used: 

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i all_ahmm_in_aims.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

#### Summary Plots
##### Across all chromosomes
Split chromosomes into five chunks each and marked any windows with at least one variant with posterior > 0.95 of homozygous _E. cordata_ ancestry.

```R
library(ggplot2)
library(RColorBrewer)

source("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/heatmaps_fun.r")

post_file_loc <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/posteriors"
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/heatmaps"
glob_mr_samples <- c("WA01", "WA03", "WA04", "WB02", "WB03", "WB04", "WC02", "WC03", "WC05", "WD04", "WE02", "WE03", "WE04", "WE05", "WF01", "WG03", "WG04", "WG05", "WH03", "WH04")
chromosomes <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

setwd(wdir)
coarsewin_filename <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/heatmaps/coarse_windows.csv"
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
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry - AIMs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
hom_coarse

het_coarse <- ggplot(coarse_mat_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry - AIMs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
het_coarse
```

![Presence of homozygous _E. cordata_ ancestry > 0.95 posterior per chromosome and sample; chromosomes 7,8 have regions shared between a large number of individuals.](heatmaps/coarse_hom_windows.png "Homozygous _E. cordata_ ancestry heatmap, AIMs only")

Similar to results with all variants, but with some former positive hits missing.
Blocks of interest for homozygous _E. cordata_ ancestry at the following based on sharing between individuals:

| Chromosome | Start (bp)  | End (bp)    | Number of Inds Sharing |
| ---------- | ----------- | ----------- | ---------------------- |
| Chr01      | 16,887,823  | 25,331,733  | 4 (16 with het)        |
| Chr07      | 1           | 10,850,526  | 8                      |
| Chr08      | 14,042,923  | 28,085,844  | 6 (14 with het)        |

##### Finer Windows within Chromosomes

```R
# Chromosome 1, Block 3
chr01b3_samples <- c("WC05", "WB02", "WD04", "WH04", "WG03", "WF01", "WG04", "WC03", "WE02", "WC02", "WB03", "WE04", "WA04", "WH03", "WE03", "WG04")
chr01b3_filename_rel <- "chr01b3.csv"
chr01b3_filename <- paste(wdir, chr01b3_filename_rel, sep = "/")
chr01b3_tab <- read.csv(chr01b3_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr01b3_df <- populate_mat(chr01b3_samples, chr01b3_tab, 0.95, post_file_loc)

chr01b3_plots <- plot_heatmap(chr01b3_df, name_table$RAPiD_ID, name_table$Accession, "Chr01", 0.95, "- Chr01 B3 - AIMs")
chr01b3_plots$het_plot

# Chromosome 1, Block 3-89
chr01b3_8_filename_rel <- "chr01b3-8.csv"
chr01b3_8_filename <- paste(wdir, chr01b3_8_filename_rel, sep = "/")
chr01b3_8_tab <- read.csv(chr01b3_8_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr01b3_8_df <- populate_mat(chr01b3_samples, chr01b3_8_tab, 0.95, post_file_loc)

chr01b3_8_plots <- plot_heatmap(chr01b3_8_df, name_table$RAPiD_ID, name_table$Accession, "Chr01", 0.95, "- Chr01 B3.8-9 - AIMs")
chr01b3_8_plots$het_plot + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))

# Chromosome 7, Block 1
chr07b1_samples <- c("WH04", "WG03", "WC03", "WE02", "WC02", "WA01", "WA04", "WH03")
chr07b1_filename_rel <- "chr07b1.csv"
chr07b1_filename <- paste(wdir, chr07b1_filename_rel, sep = "/")
chr07b1_tab <- read.csv(chr07b1_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr07b1_df <- populate_mat(chr07b1_samples, chr07b1_tab, 0.95, post_file_loc)

chr07b1_plots <- plot_heatmap(chr07b1_df, name_table$RAPiD_ID, name_table$Accession, "Chr07", 0.95, "- Chr07 B1 - AIMs")
chr07b1_plots$het_plot

# Chromosome 8, Block 2
chr08b2_samples <- c("WB02", "WD04", "WH04", "WG03", "WG05", "WE02", "WC02", "WB04", "WB03", "WE04", "WA01", "WA04", "WE03", "WE05")
chr08b2_filename_rel <- "chr08b2.csv"
chr08b2_filename <- paste(wdir, chr08b2_filename_rel, sep = "/")
chr08b2_tab <- read.csv(chr08b2_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b2_df <- populate_mat(chr08b2_samples, chr08b2_tab, 0.95, post_file_loc)

chr08b2_plots <- plot_heatmap(chr08b2_df, name_table$RAPiD_ID, name_table$Accession, "Chr08", 0.95, "- Chr08 B2 - AIMs")
chr08b2_plots$het_plot

# Chromosome 8, Block 2-7
chr08b2_7_filename_rel <- "chr08b2-7.csv"
chr08b2_7_filename <- paste(wdir, chr08b2_7_filename_rel, sep = "/")
chr08b2_7_tab <- read.csv(chr08b2_7_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))
chr08b2_7_df <- populate_mat(chr08b2_samples, chr08b2_7_tab, 0.95, post_file_loc)

chr08b2_7_plots <- plot_heatmap(chr08b2_7_df, name_table$RAPiD_ID, name_table$Accession, "Chr08", 0.95, "- Chr08 B2 - AIMs")
chr08b2_7_plots$het_plot + theme(axis.text.x = element_text(angle = 60, vjust = 0.95, hjust=1))
```

Chromosome 7 was able to be narrowed down to the blocks of interest immediately.

| Chromosome   | Start         | End           | Number of Inds Sharing |
| ------------ | ------------- | ------------- | ---------------------- |
| Chromosome 1 | 22,200,000 bp | 22,240,000 bp | 3                      |
| Chromosome 7 | 430,000 bp    | 445,000 bp    | 8                      |
| Chromosome 8 | 22,470,000 bp | 22,485,000 bp | 3 (11 with het)        |

![Presence of _E. cordata_ ancestry in Chromosome 1 Block 3; heterozygotes were scored as half in presence/absence matrix.](chr01b3-8.png "_E. cordata_ ancestry heatmap for Chromosome 1 Block 3_8")

![Presence of _E. cordata_ ancestry in Chromosome 7 Block 1; heterozygotes were scored as half in presence/absence matrix.](chr07b1.png "_E. cordata_ ancestry heatmap for Chromosome 7 Block 1")

![Presence of _E. cordata_ ancestry in Chromosome 8 Block 2; heterozygotes were scored as half in presence/absence matrix.](chr08b2-7.png "_E. cordata_ ancestry heatmap for Chromosome 8 Block 2_7")

The only region with a lot of individuals is Chromosome 7 Block 1, which is somewhat concerning given that this region is also much sparser in variant calling density.

## Low Coverage Analysis
### Generate Input files
```bash
# Run on UFRC queue system; see get_ahmm_in_aims_lc.job for more details.
# Resources used:

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/02.aims/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py "$LIST_DIR"/02.aims/ref_files.txt "$LIST_DIR"/02.aims/low_cov/coverage.txt "$LIST_DIR"/02.aims/sample_files.txt gen_dists.tab all_ahmm_in_aims_lc.tab
```

```bash
# Run on UFRC queue system; see ancestryhmm_aims_lc.job for more details.
# Resources used:

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i all_ahmm_in_aims_lc.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

### Plot results
```R
library(ggplot2)
library(RColorBrewer)

source("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/heatmaps_fun.r")

post_file_loc <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/low_cov/posteriors"
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ancestry_hmm/02.aims/low_cov/heatmaps"
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
              ggtitle("At least one locus with 95% homozygous E. cordata ancestry - AIMs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
hom_coarse

het_coarse <- ggplot(coarse_mat_df, aes(window, acc, fill = het_mat)) + 
              geom_tile(color = "#005267") + 
              scale_fill_steps(high = "#ffd966", low = "#007c9b") + 
              guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
              ggtitle("At least one locus with 95% E. cordata ancestry - AIMs") + 
              xlab("Chromosome Window") + ylab("E. globulus Sample") + 
              scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes) +
              theme(text=element_text(size=16))
het_coarse
```

