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
declare -a CHRLIST=(Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
cat Chr01/Chr01_all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt > all_fil_biallelic_globMR_all.PREDICT.BSCORRECTED.txt

for NAME in "${CHRLIST[@]}"
do
    tail -n +2 "$NAME"/"$NAME"_all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt >> all_fil_biallelic_globMR_all.PREDICT.BSCORRECTED.txt
done
```

Visualized recombination over windows and calculated genome-wide average.

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
    chr_tab <- read.table(paste(chr, "all_fil_biallelic_globMR.PREDICT.BSCORRECTED.txt", sep = "_"), header = TRUE)
    avgs <- get_avgs(chr_tab)
    winsize <- c(winsize, avgs[["winsize"]])
    nsites <- c(nsites, avgs[["nsites"]])
    recomb <- c(recomb, avgs[["recomb"]])
    graph_rwin(chr_tab, paste(chr, "_recomb_glob.png", sep = ""))
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
    graph_rwin(chr_tab, paste(chr, "_recomb_cord.png", sep = ""))
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

