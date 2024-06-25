# Local Ancestry Inference with `ELAI`
Used the local ancestry inference software [`Efficient Local Ancestry Inference`](https://github.com/haplotype/elai) to infer admixture breakpoints along genome. This software uses a two-level HMM model based on the multispecies coalescent, with the upper layer corresponding to different source populations and the lower level corresponding to haplotype variation within a population.

## Input File Preparation
Used `PLINK` to convert VCF files to `BIMBAM` format genotype files using the [`recode`](https://www.cog-genomics.org/plink/1.9/data#recode) function.

```bash
# Run on UFRC queue system; see make_bimbam.job for more details.
# Resources used: 1 Gb, 7 min

module load plink/1.90b3.39 

VCF="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"

# MR cordata
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Ecordata.fam --out "$WDIR"/cord --recode bimbam

# ref globulus
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Eglobulus_ref.fam --out "$WDIR"/glob_pure --recode bimbam

# MR globulus
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Eglobulus_MR.fam --out "$WDIR"/glob_mr --recode bimbam
```

Split the genotype files by chromosome (as per the recommendation of the `ELAI` manual).

```bash
# Run on UFRC queue system; see split_bimbam.job for more details.
# Resources used: 17 Mb, 2 sec

module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
python "$SCRIPTS_DIR"/split_bimbam.py
```

Replaced the number of SNPs in headings

```bash
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai/geno_by_chr"
cd "$WDIR"

declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

wc -l pos_Chr*.txt
declare -a SIZELIST=(787645 839980 736729 620251 620838 1079716 661655 1037845 688961 780227 774355)

rename _geno.txt _geno.tmp *

for INDEX in {0..10};
do
    head -n 1 cord_"${CHRLIST[$INDEX]}"_geno.tmp > cord_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> cord_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 cord_"${CHRLIST[$INDEX]}"_geno.tmp >> cord_"${CHRLIST[$INDEX]}"_geno.txt

    head -n 1 glob_ref_"${CHRLIST[$INDEX]}"_geno.tmp > glob_ref_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> glob_ref_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 glob_ref_"${CHRLIST[$INDEX]}"_geno.tmp >> glob_ref_"${CHRLIST[$INDEX]}"_geno.txt

    head -n 1 glob_mr_"${CHRLIST[$INDEX]}"_geno.tmp > glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 glob_mr_"${CHRLIST[$INDEX]}"_geno.tmp >> glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
done

rm *.tmp
```

## Run `ELAI`
### MAF = 0.00
Ran ELAI for all SNPs and estimated generation time since admixture = 800 (40,000 years ago / 50 years per gen). 10 separate runs were performed and results were averaged across runs as per recommendations of the developer; example given below.

```bash
# Run on UFRC queue system; see elai_maf00_g800_r1.job for more details.
# Resources used: 

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for RUN in {1..10}
do
    for CHR in "${CHRLIST[@]}"
    do
        "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r"$RUN" -C 2 -c 10 -mg 800 --exclude-miss1
    done
done
```

Averaged dosages over the 10 runs for MAF = 0.00 and generation time = 800.
```bash 
# performed on UFRC queue system; see avg_elai_maf00_g800.job for more details.
# Resources used:

module load R/4.2
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
OUTPUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai/maf00/output"

for CHR in "${CHRLIST[@]}"
do
    ls "$OUTPUT_DIR"/"$CHR"_*.ps21.txt > "$CHR"_infile_list.txt
    Rscript "$SCRIPT_DIR"/avg_elai_doses.r "$CHR"_infile_list.txt "$OUTPUT_DIR"/avg/"$CHR"_avg.ps21.txt
done
```

### MAF = 0.05
Ran ELAI for MAF = 0.05 and estimated generation time since admixture = 800.

```bash
# Run on UFRC queue system; see elai_maf05_g800.job for more details.
# Resources used: 3 Gb, 4 days

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for RUN in {1..10}
do
    for CHR in "${CHRLIST[@]}"
    do
        "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r"$RUN" -C 2 -c 10 -mg 800 -exclude-maf 0.05
    done
done
```

Averaged dosages over the 10 runs for MAF = 0.00 and generation time = 800.
```bash 
# performed on UFRC queue system; see avg_elai_maf05_g800.job for more details.
# Resources used: 5 Gb, 6 days

module load R/4.2
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
OUTPUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai/maf05/output"

for CHR in "${CHRLIST[@]}"
do
    ls "$OUTPUT_DIR"/"$CHR"_*.ps21.txt > file_lists/"$CHR"_infile_list.txt
    Rscript "$SCRIPT_DIR"/avg_elai_doses.r file_lists/"$CHR"_infile_list.txt "$OUTPUT_DIR"/avg/"$CHR"_avg.ps21.txt
done
```

Plotted heatmaps of average _E. cordata_ dosage in each putative introgressed _E. globulus_ individual using custom `R` functions (see `heatmaps_elai_fun.r`).

```R
library(ggplot2)
library(RColorBrewer)

source("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/elai/heatmaps_elai_fun.r")

infile_loc <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/elai/maf05/results"
wdir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/elai/maf05/heatmaps"

sample_order <- read.table(paste(wdir, "..", "..", "Eglobulus_MR_inds.txt", sep = "/"))$V1
chromosomes <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

setwd(wdir)
coarsewin_filename <- paste(wdir, "..", "..", "coarse_windows.csv", sep = "/")
coarsewin_tab <- read.csv(coarsewin_filename, header = FALSE, col.names = c("chrom", "start", "end"), colClasses = c("character", "integer", "integer"))

# populate matrix
coarse_mat_df <- populate_elai_mat_val(infile_loc, sample_order, coarsewin_tab)

# plot heatmaps
x_axis_breaks <- c("Chr01_16887823_25331733", "Chr02_20331353_30497028", "Chr03_26218899_39328347", "Chr04_15439735_23159601"
, "Chr05_25167991_37751985", "Chr06_20856265_31284396", "Chr07_21701053_32551578", "Chr08_28085845_42128766", "Chr09_15320131_22980195", "Chr10_15489065_23233596", "Chr11_16822585_25233876")

name_table_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
name_table <- read.csv(name_table_name, header = TRUE, as.is = TRUE)
label_order <- match(coarse_mat_df$sample, name_table$RAPiD_ID)
coarse_mat_df$acc <- name_table$Accession[label_order]

dose_coarse <- ggplot(coarse_mat_df, aes(window, acc, fill = dosage)) + 
               geom_tile(color = "#000000") + 
               scale_fill_steps(high = "#000000", low = "#ffffff") + 
               guides(fill = guide_coloursteps(title = NULL, show.limits = TRUE)) +
               ggtitle("Maximum dosage of E. cordata ancestry per window") + 
               xlab("Chromosome Window") + ylab("E. globulus Sample") + 
               scale_x_discrete(breaks = x_axis_breaks, labels = chromosomes)
png("genwide_maf05_display.png", width = 5000, height = 1750)
print(dose_coarse +
      theme(axis.text = element_text(size = 70),
            axis.title = element_text(size = 78),
            plot.title = element_text(size = 95),
            legend.text = element_text(size = 45),
            legend.title = element_text(size = 50),
            legend.key.size = unit(3, "cm"))
      )
dev.off()
```
