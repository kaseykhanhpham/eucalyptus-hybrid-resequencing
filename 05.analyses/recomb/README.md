# Recombination Landscape
Used [`ReLERNN`](https://github.com/kr-colab/ReLERNN), an deep learning ML approach, to estimate recombination across the genome for _E. globulus_ and _E. cordata_ at Meehan Range.

Filtered VCF to biallelic SNPs for one sample group.
```bash
module load bcftools/1.15
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
cd "$VCF_DIR"

bcftools view -O v -o all_fil_biallelic_globMR.vcf --samples-file /blue/soltis/kasey.pham/euc_hyb_reseq/analyses/recomb/glob/Eglobulus_MR_samples.txt --threads 16 all_fil_biallelic.vcf.gz

bcftools view -O v -o all_fil_biallelic_cord.vcf --samples-file /blue/soltis/kasey.pham/euc_hyb_reseq/analyses/recomb/cord/Ecordata_samples.txt --threads 16 all_fil_biallelic.vcf.gz
```
## Estimating _E. globulus_ recombination
Created training set using `ReLERNN_SIMULATE`, sticking to defaults for the number of sets made for each step of training. The manual of `ReLERNN` also recommended not using information from a demographic plot unless very confident in the results, so I abstained from doing that.

```bash
# Run in UFRC's job queue system; see relernn_simulate_glob.job for more details.
# Resources used: 7 Gb, 1 min

module load relernn/1.0.0
module load cuda/12.2.2
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

# mutation rate and generation time from Silva-Junior and Grattapaglia 2015 New Phytol.

for NAME in "${VCFLIST[@]}"
do
    ReLERNN_SIMULATE -v "$VCF_DIR"/all_fil_biallelic_globMR.vcf -g ../"$NAME".bed -d "$NAME" -u "4.93e-8" -l 10 -t 8 --unphased
done
```

Trained recurrent network using simulated datasets. Used defaults for training periods.
```bash
# Run in UFRC's job queue system; see relernn_train_glob.job for more details.
# Resources used: 76 Gb, 17 hrs

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
# Resources used: 54 Gb, 2 hrs

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

Visualized recombination over windows and calculated genome-wide average.

```R
source("../recomb_graph_funs.r")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

# create storage vectors
winsize <- c()
nsites <- c()
recomb <- c()

# loop through each chromosome-specific table and average the window size, number of sites per window, and recombination rate.
for(chr in chr_list){
    chr_tab <- read.table(paste("all_fil_biallelic_globMR_", chr, ".PREDICT.BSCORRECTED.txt", sep = ""), header = TRUE)
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
write.table(out_tab, "recomb_avgs_glob.tab", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
```

| Chromosome  | Window Size | Num Sites | Recombination (/bp*gen) |
| ----------- | ----------- | --------- | ----------------------- |
| Chr01       | 301492.88   | 729.29    | 1.717e-08               |
| Chr02       | 276153.82   | 582.07    | 1.965e-08               |
| Chr03       | 212000      | 472.61    | 1.641e-08               |
| Chr04       | 364093.33   | 764.36    | 1.513e-08               |
| Chr05       | 329000      | 610.26    | 2.222e-08               |
| Chr06       | 152892.56   | 384.37    | 1.657e-08               |
| Chr07       | 268452.74   | 548.43    | 2.271e-08               |
| Chr08       | 378000      | 844.02    | 1.668e-08               |
| Chr09       | 273530.11   | 670.81    | 1.597e-08               |
| Chr10       | 304000      | 762.04    | 2.198e-08               |
| Chr11       | 288000      | 645.63    | 2.014e-08               |
| genome-wide | 286146.86   | 637.62    | 1.860e-08               |


## Estimate _E. cordata_ recombination
Repeated the above steps for _E. cordata_-only VCF file. Then renamed and re-organized output tabs for processing.

```bash
declare -a VCFLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)
for NAME in "${VCFLIST[@]}"; do rename all_fil_biallelic_cord all_fil_biallelic_cord_"$NAME" "$NAME"/*; done
```

Calculated chromosome-wide and genome-wide average recombination and plotted recombination windows.
```bash
source("../recomb_graph_funs.r")
chr_list <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")

# create storage vectors
winsize <- c()
nsites <- c()
recomb <- c()

# loop through each chromosome-specific table and average the window size, number of sites per window, and recombination rate.
for(chr in chr_list){
    chr_tab <- read.table(paste("all_fil_biallelic_cord_", chr, ".PREDICT.BSCORRECTED.txt", sep = ""), header = TRUE)
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
| Chr01       | 301492.88   | 729.29    | 2.903e-08               |
| Chr02       | 276153.82   | 582.07    | 2.855e-08               |
| Chr03       | 212000      | 472.61    | 1.756e-08               |
| Chr04       | 364093.33   | 764.36    | 3.412e-08               |
| Chr05       | 329000      | 610.26    | 8.905e-09               |
| Chr06       | 152892.56   | 384.37    | 2.293e-08               |
| Chr07       | 268452.74   | 548.43    | 1.894e-08               |
| Chr08       | 378000      | 844.02    | 2.802e-08               |
| Chr09       | 273530.11   | 670.81    | 2.229e-08               |
| Chr10       | 304000      | 762.04    | 2.380e-08               |
| Chr11       | 288000      | 645.63    | 1.967e-08               |
| genome-wide | 286146.86   | 637.62    | 2.307e-08               |
