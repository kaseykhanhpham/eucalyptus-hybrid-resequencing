# Recombination Landscape
Used [`ReLERNN`](https://github.com/kr-colab/ReLERNN), an deep learning ML approach, to estimate recombination across the genome for _E. globulus_ and _E. cordata_ at Meehan Range.

Filtered VCF to biallelic SNPs for one sample group.
```bash
module load bcftools/1.15
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
cd "$VCF_DIR"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

bcftools view -O v -o all_fil_biallelic_globMR.vcf --samples-file "$LIST_DIR"/Eglobulus_MR.txt --threads 16 all_fil_biallelic.vcf.gz

bcftools view -O v -o all_fil_biallelic_cord.vcf --samples-file "$LIST_DIR"/Ecordata.txt --threads 16 all_fil_biallelic.vcf.gz
```
## Estimating _E. globulus_ recombination
### ReLERNN
Created training set using `ReLERNN_SIMULATE`, sticking to defaults for the number of sets made for each step of training. The manual of `ReLERNN` also recommended not using information from a demographic plot unless very confident in the results, so I abstained from doing that.

```bash
# Run in UFRC's job queue system; see relernn_simulate_glob.job for more details.
# Resources used: 

module load relernn/1.0.0
module load cuda/12.2.2
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

# mutation rate (4.93e-8 per base pair per generation) from Silva-Junior and Grattapaglia 2015 New Phytol. generation time (50 yrs per generation) from discussion with Potts/Vaillancourt lab.

for NAME in "${VCFLIST[@]}"
do
    ReLERNN_SIMULATE -v "$VCF_DIR"/all_fil_biallelic_globMR.vcf -g ../"$NAME".bed -d "$NAME" -u "4.93e-8" -l 50 -t 8 --unphased
done
```

Trained recurrent network using simulated datasets. Used defaults for training periods.
```bash
# Run in UFRC's job queue system; see relernn_train_glob.job for more details.
# Resources used: 100 Gb, 1 day

module load relernn/1.0.0
module load cuda/12.2.2
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/apps/compilers/cuda/12.2.2

for NAME in "${VCFLIST[@]}"
do
    ReLERNN_TRAIN -d "$NAME" -t 8
done
```

Predicted recombination in windows of 40 SNPs across genome using trained model.
```bash
# Run in UFRC's job queue system; see relernn_predict_glob.job for more details.
# Resources used: 780 Mb, 2 min

module load relernn/1.0.0
module load cuda/12.2.2
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/apps/compilers/cuda/12.2.2

for NAME in "${VCFLIST[@]}"
do
    ReLERNN_PREDICT -v "$VCF_DIR"/all_fil_biallelic_globMR.vcf -d "$NAME" --unphased --minSites 40
done
```

Estimated confidence intervals using simulated replicates of each recombination window.
```bash
# Run in UFRC's job queue system; see relernn_bscorrect_glob.job for more details.
# Resources used: 53 Gb, 3 hrs

module load relernn/1.0.0
module load cuda/12.2.2
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/apps/compilers/cuda/12.2.2

for NAME in "${VCFLIST[@]}"
do
    ReLERNN_BSCORRECT -d "$NAME" -t 8
done
```

Renamed ReLERNN outputs by chromosome.
```bash
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${CHRLIST[@]}"
do
    rename all_fil_biallelic_globMR "$NAME"_all_fil_biallelic_globMR "$NAME"/*
done
```

Consolidated window files into one.
```bash
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
head -n 1 Chr01/Chr01_all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt > all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt

for NAME in "${CHRLIST[@]}"
do
    tail -n +2 "$NAME"/"$NAME"_all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt >> all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt
done
```
### Plotting
Visualized recombination over windows and calculated genome-wide average.

```R
source("../../recomb_graph_funs.r")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

# create storage vectors
winsize <- c()
nsites <- c()
recomb <- c()

# loop through each chromosome-specific table and average the window size, number of sites per window, and recombination rate.
for(chr in chr_list){
    chr_tab <- read.table(paste(chr, "all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt", sep = "_"), header = TRUE)
    avgs <- get_avgs(chr_tab)
    winsize <- c(winsize, avgs[["winsize"]])
    nsites <- c(nsites, avgs[["nsites"]])
    recomb <- c(recomb, avgs[["recomb"]])
    graph_rwin_glob(chr_tab, paste(chr, "_recomb_glob.png", sep = ""), 0, 6.5e-08)
}

# append genome-wide averages
chr_list <- c(chr_list, "genome-wide")
avg_winsize <- mean(winsize)
winsize <- c(winsize, avg_winsize)
avg_nsites <- mean(nsites)
nsites <- c(nsites, avg_nsites)
avg_recomb <- mean(recomb)
recomb <- c(recomb, avg_recomb)

out_tab <- data.frame(chr_list, winsize, nsites, recomb)
write.table(out_tab, "../recomb_avgs_glob.tab", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```

| Chromosome  | Window Size | Num Sites | Recombination (/bp*gen) |
| ----------- | ----------- | --------- | ----------------------- |
| Chr01       | 154624      | 614       | 1.138e-08               |
| Chr02       | 179544      | 645       | 1.731e-08               |
| Chr03       | 120000      | 482       | 1.440e-08               |
| Chr04       | 177863      | 637       | 1.832e-08               |
| Chr05       | 169000      | 581       | 2.180e-08               |
| Chr06       | 76000       | 323       | 1.642e-08               |
| Chr07       | 113000      | 428       | 1.330e-08               |
| Chr08       | 168000      | 659       | 1.646e-08               |
| Chr09       | 127614      | 529       | 1.886e-08               |
| Chr10       | 153000      | 641       | 1.679e-08               |
| Chr11       | 164000      | 623       | 2.105e-08               |
| genome-wide | 145695      | 560       | 1.692e-08               |

(Window size and number of sites rounded to the nearest whole number.)

## Estimate _E. cordata_ recombination
Repeated the above steps for _E. cordata_-only VCF file. Then renamed and re-organized output tabs for processing.

```bash
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
for NAME in "${VCFLIST[@]}"; do rename all_fil_biallelic_cord "$NAME"_all_fil_biallelic_cord "$NAME"/*; done
```

### Plotting
Calculated chromosome-wide and genome-wide average recombination and plotted recombination windows.
```R
library(ggplot2)
source("../../recomb_graph_funs.r")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

# create storage vectors
winsize <- c()
nsites <- c()
recomb <- c()

# loop through each chromosome-specific table and average the window size, number of sites per window, and recombination rate.
for(chr in chr_list){
  chr_tab <- read.table(paste(chr, "all_fil_biallelic_cord.PREDICT.BSCORRECTED.txt", sep = "_"), header = TRUE)
  avgs <- get_avgs(chr_tab)
  winsize <- c(winsize, avgs[["winsize"]])
  nsites <- c(nsites, avgs[["nsites"]])
  recomb <- c(recomb, avgs[["recomb"]])
  graph_rwin_cord(chr_tab, paste(chr, "_recomb_cord.png", sep = ""), 0, 6.5e-08)
}

# append genome-wide averages
chr_list <- c(chr_list, "genome-wide")
avg_winsize <- mean(winsize)
winsize <- c(winsize, avg_winsize)
avg_nsites <- mean(nsites)
nsites <- c(nsites, avg_nsites)
avg_recomb <- mean(recomb)
recomb <- c(recomb, avg_recomb)

out_tab <- data.frame(chr_list, winsize, nsites, recomb)
write.table(out_tab, "recomb_avgs_cord.tab", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```

| Chromosome  | Window Size | Num Sites | Recombination (/bp*gen) |
| ----------- | ----------- | --------- | ----------------------- |
| Chr01       | 154624      | 614       | 2.002e-08               |
| Chr02       | 179544      | 645       | 2.182e-08               |
| Chr03       | 120000      | 482       | 1.671e-08               |
| Chr04       | 177863      | 637       | 1.833e-08               |
| Chr05       | 169000      | 581       | 9.520e-09               |
| Chr06       | 76000       | 323       | 1.137e-08               |
| Chr07       | 113000      | 428       | 1.165e-08               |
| Chr08       | 168000      | 659       | 2.441e-08               |
| Chr09       | 127614      | 529       | 1.553e-08               |
| Chr10       | 153000      | 641       | 2.357e-08               |
| Chr11       | 164000      | 623       | 2.613e-08               |
| genome-wide | 145695      | 560       | 1.810e-08               |

Consolidated window files into one.
```bash
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
head -n 1 Chr01/Chr01_all_fil_biallelic_cord.PREDICT.BSCORRECTED.txt > all_fil_biallelic_cord_all.PREDICT.BSCORRECTED.txt

for NAME in "${CHRLIST[@]}"
do
    tail -n +2 "$NAME"/"$NAME"_all_fil_biallelic_cord.PREDICT.BSCORRECTED.txt >> all_fil_biallelic_cord_all.PREDICT.BSCORRECTED.txt
done
```

## Recombination versus LD and gene density
Tested correlation of recombination windows in both taxa to gene density and sliding window LD.
```R
# CORRELATION BETWEEN LD AND RECOMBINATION
# define functions for processing tables
get_ld_pos <- function(file_name_list) {
  ld_interv <- sapply(strsplit(file_name_list, "_"), function (x) x[2])
  ld_starts <- as.numeric(sapply(strsplit(ld_interv, "-"), function (x) x[1]))
  ld_ends <- as.numeric(sapply(strsplit(ld_interv, "-"), function (x) x[2]))
  return(list(starts = ld_starts, ends = ld_ends))
}

avg_ld_recwins <- function(rec_tab, ld_tab) {
  # extract LD window positions
  ld_pos <- get_ld_pos(ld_tab$file_name)
  # get LD windows overlapping recombination windows and average LD from them
  ld_in_rec_wins <- c()
  for (i in c(1:nrow(rec_tab))) {
    totally_within <- intersect(which(ld_pos[["starts"]] > rec_tab[i, "start"]),
                                which(ld_pos[["ends"]] < rec_tab[i, "end"]))
    start_overl <- intersect(which(ld_pos[["starts"]] < rec_tab[i, "start"]),
                             which(ld_pos[["ends"]] > rec_tab[i, "start"]))
    end_overl <- intersect(which(ld_pos[["starts"]] < rec_tab[i, "end"]),
                           which(ld_pos[["ends"]] > rec_tab[i, "end"]))
    all_overl <- sort(unique(c(totally_within, start_overl, end_overl)))
    avg_ld <- mean(ld_tab[all_overl, "ld"], na.rm = TRUE)
    ld_in_rec_wins <- c(ld_in_rec_wins, avg_ld)
  }
  # return vector of average LD corresponding to recombination table rows
  return(ld_in_rec_wins)
}

# E. globulus
glob_ld <- read.table("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ld/glob/curve_fit/windows/glob_all_windows_ld.txt", header = TRUE)
glob_rec <- read.table("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/recomb/glob/all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt", header = TRUE)

# calculate average LD overlapping with recomb windows
glob_ld_recwins <- avg_ld_recwins(glob_rec, glob_ld)

# plot LD vs. RECOMBINATION RATE
ld_rec_lm <- lm(glob_ld_recwins ~ glob_rec$recombRate) # P = 0.000774, r2 = 0.003174
plot(glob_rec$recombRate, glob_ld_recwins, pch = 16, col = "#13BDD7",
     main = "E. globulus LD vs. Recombination Rate",
     xlab = "Recombination Rate (1/bp*gen)", ylab = "LD (bp)", cex.axis = 0.9)
abline(ld_rec_lm, lwd = 3, col = "#006185")
text(x = 4e-08, y = 12500, "P = 0.000774", adj = 0)
text(x = 4e-08, y = 11800, "r^2 = 0.003174", adj = 0)

# E. cordata
cord_ld <- read.table("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/ld/cord/curve_fit/windows/cord_all_windows_ld.txt", header = TRUE)
cord_rec <- read.table("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/recomb/cord/all_fil_biallelic_cord.PREDICT.BSCORRECTED.txt", header = TRUE)

# calculate average LD overlapping with recomb windows
cord_ld_recwins <- avg_ld_recwins(cord_rec, cord_ld)

# plot LD vs. RECOMBINATION RATE
ld_rec_lm <- lm(cord_ld_recwins ~ cord_rec$recombRate) # P = 0.257, r2 = 0.00008779
plot(cord_rec$recombRate, cord_ld_recwins, pch = 16, col = "#fdc51e",
     main = "E. cordata LD vs. Recombination Rate",
     xlab = "Recombination Rate (1/bp*gen)", ylab = "LD (bp)", cex.axis = 0.9)
abline(ld_rec_lm, lwd = 3, col = "#d78700")
text(x = 4e-08, y = 29000, "P = 0.257", adj = 0)
text(x = 4e-08, y = 27500, "r^2 = 0.00008779", adj = 0)

# Get gene density in recombination windows
annot <- read.table("C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/references/Eglobulus_genome_X46/EGLOB-X46.v1.0.annotation.gff", header = FALSE, sep = "\t")
annot <- annot[which(annot$V3 == "gene"),]
# define function to iterate through recombination table and extract genes overlapping each window
get_gene_dens <- function(annot, rec_tab) {
  gene_dens <- c()
  for (i in c(1:nrow(rec_tab))) {
    totally_within <- intersect(which(annot$V4 > rec_tab[i, "start"]),
                                which(annot$V5 < rec_tab[i, "end"]))
    start_overl <- intersect(which(annot$V4 < rec_tab[i, "start"]),
                             which(annot$V5 > rec_tab[i, "start"]))
    end_overl <- intersect(which(annot$V4 < rec_tab[i, "end"]),
                           which(annot$V5 > rec_tab[i, "end"]))
    all_overl <- sort(unique(c(totally_within, start_overl, end_overl)))
    win_dens <- length(all_overl)/(rec_tab[i, "end"] - rec_tab[i, "start"])
    gene_dens <- c(gene_dens, win_dens)
  }
  return(gene_dens)
}

# E. globulus RECOMBINATION vs. GENE DENSITY
glob_gene_dens <- get_gene_dens(annot, glob_rec)
glob_gene_rec_lm <- lm(glob_gene_dens ~ glob_rec$recombRate)
plot(glob_rec$recombRate, glob_gene_dens, pch = 16, col = "#13BDD7",
     main = "E. globulus Gene Density vs. Recombination Rate",
     xlab = "Recombination Rate (1/bp*gen)", ylab = "Gene Density (1/bp)")
abline(glob_gene_rec_lm, lwd = 3, col = "#006185")
text(x = 4e-08, y = 0.0011, "P = 0.888", adj = 0)
text(x = 4e-08, y = 0.00105, "r^2 = -0.0003024", adj = 0)
# E. globulus LD vs. gene density
glob_gene_ld_lm <- lm(glob_gene_dens ~ glob_ld_recwins)
plot(glob_ld_recwins, glob_gene_dens, pch = 16, col = "#13BDD7",
     main = "E. globulus Gene Density vs. LD",
     xlab = "LD (bp)", ylab = "Gene Density (1/bp)")
abline(glob_gene_ld_lm, lwd = 3, col = "#006185")
text(x = 10000, y = 0.0011, "P = 1.602e-07", adj = 0)
text(x = 10000, y = 0.00103, "r^2 = 0.008133", adj = 0)

# E. cordata RECOMBINATION vs. GENE DENSITY
cord_gene_dens <- get_gene_dens(annot, cord_rec)
cord_gene_rec_lm <- lm(cord_gene_dens ~ cord_rec$recombRate)
plot(cord_rec$recombRate, cord_gene_dens, pch = 16, col = "#fdc51e",
     main = "E. cordata Gene Density vs. Recombination Rate",
     xlab = "Recombination Rate (1/bp*gen)", ylab = "Gene Density (1/bp)")
abline(cord_gene_rec_lm, lwd = 3, col = "#d78700")
text(x = 4e-08, y = 0.0011, "P = 0.407", adj = 0)
text(x = 4e-08, y = 0.00105, "r^2 = -9.604e-05", adj = 0)
# E. cordata LD vs. gene density
cord_gene_ld_lm <- lm(cord_gene_dens ~ cord_ld_recwins)
plot(cord_ld_recwins, cord_gene_dens, pch = 16, col = "#fdc51e",
     main = "E. cordata Gene Density vs. LD",
     xlab = "LD (bp)", ylab = "Gene Density (1/bp)")
abline(cord_gene_ld_lm, lwd = 3, col = "#d78700")
text(x = 24700, y = 0.0011, "P =  0.9276", adj = 0)
text(x = 24700, y = 0.00105, "r^2 = -0.000306", adj = 0)
```
