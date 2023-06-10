# Sliding Window Genome Statistics

## Nucleotide Diversity and Divergence

Used [`pixy`](https://pixy.readthedocs.io/en/latest) to calculate pi (nucleotide diversity), dxy (absolute divergence), and FST (relative divergence) across the genome. First, removed the outgroup _E. grandis_ sample from VCFs and the outlier sample WF03/1051 for checking.

```bash
# Done on UFRC queue system; see remove_outgr_outl.job for more details.
# Resources used: ??? slurm fritzed out on me.

module load bcftools/1.15

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
OUTGR_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outgr"
OUTL_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outl"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    bcftools view --samples ^SRR10339635 --output-type z --output "$OUTGR_DIR"/"$NAME"_fil_no_outgr.vcf.gz --threads 12 "$VCF_DIR"/"$NAME"_fil.vcf.gz
    bcftools view --samples ^WF03,SRR10339635 --output-type z --output "$OUTL_DIR"/"$NAME"_fil_no_outl.vcf.gz --threads 12 "$VCF_DIR"/"$NAME"_fil.vcf.gz
done
```

Moved to `/orange` and indexed with `tabix`.

```bash
module load bcftools/1.15
OUTGR_DIR="/orange/soltis/kasey.pham/eucalyptus_hybrid_resequencing/06.snp_filtering/no_outgroup"
OUTL_DIR="/orange/soltis/kasey.pham/eucalyptus_hybrid_resequencing/06.snp_filtering/no_outlier"

cd "$OUTGR_DIR"
ls *.vcf.gz | while read FILE; do tabix "$FILE"; done

cd "$OUTL_DIR"
ls *.vcf.gz | while read FILE; do tabix "$FILE"; done
```

Created symlinks of all to `/blue` before proceeding.

Ran `pixy` for all samples except the _E. grandis_ outgroup.

```bash
# Done on UFRC queue system; see pixy.job for more details.
# Resources used: 1.3 Gb, 7 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outgr"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

conda activate "$ENV_DIR"/pixy

for NAME in "${VCFLIST[@]}"
do
    pixy --stats pi fst dxy --vcf "$VCF_DIR"/"$NAME"_fil_no_outgr.vcf.gz --populations "$LIST_DIR"/populations.txt --window_size 5000 --n_cores 12 --output_prefix "$NAME" --fst_type wc
done
```

Ran `pixy` for all samples except WF03/1051 and the `E. grandis` outgroup.

```bash
# Done on UFRC queue system; see pixy_outlier.job for more details.
# Resources used: 1.3 Gb, 7 min

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/no_outl"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/outlier_check"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

conda activate "$ENV_DIR"/pixy

for NAME in "${VCFLIST[@]}"
do
    pixy --stats pi fst dxy --vcf "$VCF_DIR"/"$NAME"_fil_no_outl.vcf.gz --populations "$LIST_DIR"/populations_no_outl.txt --window_size 5000 --n_cores 12 --output_prefix "$NAME" --fst_type wc
done
```

### Compare `pixy` results with and without outlier WF03/1051
Consolidate results:

```bash
# all but outgroup
OUTGR_DIR="/orange/soltis/kasey.pham/eucalyptus_hybrid_resequencing/06.snp_filtering/no_outgroup"
OUTL_DIR="/orange/soltis/kasey.pham/eucalyptus_hybrid_resequencing/06.snp_filtering/no_outlier"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

cd "$OUTGR_DIR"
touch all_pi.txt
touch all_dxy.txt
touch all_fst.txt

for NAME in "${VCFLIST[@]}"
do
    cat "$NAME"_pi.txt >> all_pi.txt
    cat "$NAME"_dxy.txt >> all_dxy.txt
    cat "$NAME"_fst.txt >> all_fst.txt
done

cd "$OUTL_DIR"
touch all_outl_pi.txt
touch all_outl_dxy.txt
touch all_outl_fst.txt

for NAME in "${VCFLIST[@]}"
do
    cat "$NAME"_pi.txt >> all_outl_pi.txt
    cat "$NAME"_dxy.txt >> all_outl_dxy.txt
    cat "$NAME"_fst.txt >> all_outl_fst.txt
done
```

Analyzed results in `R` on local computer.
```R
library(ggplot2)
library(dplyr)

working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/genome_scan/pixy"
setwd(working_dir)

# import pixy files
all_pi <- read.table("all_pi.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)
all_dxy <- read.table("all_dxy.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)
all_fst <- read.table("all_fst.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)

outl_pi <- read.table("all_outl_pi.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)
outl_dxy <- read.table("all_outl_dxy.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)
outl_fst <- read.table("all_outl_fst.txt", header = TRUE, sep = "\t", na.strings="NA", as.is = TRUE)

# cast column classes
pi_to_cast <- c("avg_pi", "no_sites", "count_diffs", "count_comparisons", "count_missing")
dxy_to_cast <- c("avg_dxy", "no_sites", "count_diffs", "count_comparisons", "count_missing")
fst_to_cast <- c("avg_wc_fst", "no_snps")

for(col_name in pi_to_cast){
    all_pi[,col_name] <- as.numeric(all_pi[,col_name])
    outl_pi[,col_name] <- as.numeric(outl_pi[,col_name])
}

for(col_name in dxy_to_cast){
    all_dxy[,col_name] <- as.numeric(all_dxy[,col_name])
    outl_dxy[,col_name] <- as.numeric(outl_dxy[,col_name])
}

for(col_name in fst_to_cast){
    all_fst[,col_name] <- as.numeric(all_fst[,col_name])
    outl_fst[,col_name] <- as.numeric(outl_fst[,col_name])
}

# filter by number of sites used to calculate in window
all_pi_fil <- all_pi[which(all_pi$no_sites > 40),]
outl_pi_fil <- outl_pi[which(outl_pi$no_sites > 40),]

all_dxy_fil <- all_dxy[which(all_dxy$no_sites > 40),]
outl_dxy_fil <- outl_dxy[which(outl_dxy$no_sites > 40),]

all_fst_fil <- all_fst[which(all_fst$no_snps > 15),]
outl_fst_fil <- outl_fst[which(outl_fst$no_snps > 15),]

# compare genome-wide values of pi and dxy for datasets with and without the outlier, WF03/1051

# PI
allpops_pi <- sum(all_pi_fil$count_diffs)/sum(all_pi_fil$count_comparisons) # 0.02827
allpops_pi_outl <- sum(outl_pi_fil$count_diffs)/sum(outl_pi_fil$count_comparisons) # 0.02843

glob_pi <- sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_comparisons"]) # 0.02790
glob_pi_outl <- sum(outl_pi_fil[which(outl_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_diffs"])/sum(outl_pi_fil[which(outl_pi_fil$pop %in% c("glob_MR", "glob_pure")),"count_comparisons"]) # 0.02809

ref_pi <- sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(all_pi_fil[which(all_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02802
ref_pi_outl <- sum(outl_pi_fil[which(outl_pi_fil$pop == "glob_pure"),"count_diffs"])/sum(outl_pi_fil[which(outl_pi_fil$pop == "glob_pure"),"count_comparisons"]) # 0.02938

# DXY
allpops_dxy <- sum(all_dxy_fil$count_diffs)/sum(all_dxy_fil$count_comparisons) # 0.03554
allpops_dxy_outl <- sum(outl_dxy_fil$count_diffs)/sum(outl_dxy_fil$count_comparisons) # 0.03610

ref_mr_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_MR" & all_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.02794
ref_mr_dxy_outl <- sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_MR" & outl_dxy_fil$pop2 == "glob_pure"),"count_diffs"])/sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_MR" & outl_dxy_fil$pop2 == "glob_pure"),"count_comparisons"]) # 0.02857

ref_cord_dxy <- sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(all_dxy_fil[which(all_dxy_fil$pop1 == "glob_pure" & all_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.04022
ref_cord_dxy_outl <- sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_pure" & outl_dxy_fil$pop2 == "cord_MR"),"count_diffs"])/sum(outl_dxy_fil[which(outl_dxy_fil$pop1 == "glob_pure" & outl_dxy_fil$pop2 == "cord_MR"),"count_comparisons"]) # 0.04089
```

Whole-genome pi comparisons:
| group      | pi with WF03/1051 | pi without WF03/1051 |
| ---------- | ----------------- | -------------------- |
| all        | 0.02827           | 0.02843              |
| _E. glob_  | 0.02790           | 0.02809              |
| _glob_ ref | 0.02802           | 0.02938              |

Whole-genome dxy comparisons:
| group 1   | group 2    | dxy with WF03/1051 | dxy without WF03/1051 |
| --------- | ---------- | ------------------ | --------------------- |
| all       | all        | 0.03554            | 0.03610               |
| _glob_ MR | _glob_ ref | 0.02794            | 0.02857               |
| _cord_    | _glob_ ref | 0.04022            | 0.04089               |

Inclusion of WF03/1051 changed the estimated genome-wide pi for the reference group alone the most, as expected. Values of dXY remained mostly not too strongly different. The numbers don't seem so extremely different that I would be warranted in excluding the sample...

### Identify outlier windows


## Patterson's D and associated statistics

Filtered SNP set to biallelic SNPs only using `vcftools`.

```bash
# Run on UFRC queue system; see get_biallelic.job for more details.
# Resources used: 10 Mb, 2 min

module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > "$VCFDIR"/all_fil_biallelic.vcf.gz
```

Calculated fD and fDM in sliding windows of 40 viable SNPs at a time using [`Dsuite`](https://github.com/millanek/Dsuite).

```bash
# Run on UFRC queue system; see dsuite.job for more details.
# Resources used: 4 Mb, 8 min

module load gcc/12.2.0

DSUITE_DIR="/blue/soltis/kasey.pham/bin/Dsuite/Build"
# My install is compiled on GCC 12.2.0, versus the HPC install on 9.3.0
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"

"$DSUITE_DIR"/Dsuite Dinvestigate -w 40,20 -g "$VCF_DIR"/all_fil_biallelic.vcf.gz SETS.txt test_trios.txt
```

### Get outlier windows


## Tajima's D and Transition/Transversion
Calculated Tajima's D and Ts/Tv rate in 5000bp sliding windows using `vcftools`.

```bash
# Run on UFRC queue system; see tajd_windows.job for more details.
# Resources used: 5 Mb, 3 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Ecordata.txt --TajimaD 5000 --out cord_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_MR.txt --TajimaD 5000 --out glob_mr_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_ref.txt --TajimaD 5000 --out glob_ref_"$NAME"
done
```

```bash
# Run on UFRC queue system; see tstv_windows.job for more details.
# Resources used: 5 Mb, 2 min

module load vcftools/0.1.16

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/vcftools"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Ecordata.txt --TsTv 5000 --out cord_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_MR.txt --TsTv 5000 --out glob_mr_"$NAME"
    vcftools --gzvcf "$VCF_DIR"/"$NAME"_fil.vcf.gz --keep Eglobulus_ref.txt --TsTv 5000 --out glob_ref_"$NAME"
done
```

Merged chromosome files.

```bash

```

### Get outlier windows
