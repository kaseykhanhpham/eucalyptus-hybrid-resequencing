# Linkage Disequilibrium
## Allele Phasing
Index VCF files
```bash
module load samtools/0.1.18
ls *.bam | while read FILE
do
    tabix "$FILE"
done
```

Ran WhatsHap
```bash
# Ran in UFRC queue system; see whatshap.job for more details.
# Resources Used: 11 Gb, 2 days

module load whatshap/1.1
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
BAM_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${VCFLIST[@]}"
do
    whatshap phase -o "$NAME"_hapmarked.vcf --reference "$REF_DIR"/EGLOB-X46.v1.0.fa --tag PS "$VCF_DIR"/"$NAME"_fil.vcf.gz "$BAM_DIR"/S10_marked.bam "$BAM_DIR"/S2_marked.bam "$BAM_DIR"/S337_marked.bam "$BAM_DIR"/S356_marked.bam "$BAM_DIR"/S11_marked.bam "$BAM_DIR"/S30_marked.bam "$BAM_DIR"/S338_marked.bam "$BAM_DIR"/S357_marked.bam "$BAM_DIR"/S12_marked.bam "$BAM_DIR"/S31_marked.bam "$BAM_DIR"/S339_marked.bam "$BAM_DIR"/S358_marked.bam "$BAM_DIR"/S13_marked.bam "$BAM_DIR"/S320_marked.bam "$BAM_DIR"/S33_marked.bam "$BAM_DIR"/S359_marked.bam "$BAM_DIR"/S14_marked.bam "$BAM_DIR"/S321_marked.bam "$BAM_DIR"/S340_marked.bam "$BAM_DIR"/S35_marked.bam "$BAM_DIR"/S15_marked.bam "$BAM_DIR"/S322_marked.bam "$BAM_DIR"/S341_marked.bam "$BAM_DIR"/S36_marked.bam "$BAM_DIR"/S16_marked.bam "$BAM_DIR"/S323_marked.bam "$BAM_DIR"/S342_marked.bam "$BAM_DIR"/S37_marked.bam "$BAM_DIR"/S17_marked.bam "$BAM_DIR"/S324_marked.bam "$BAM_DIR"/S343_marked.bam "$BAM_DIR"/S38_marked.bam "$BAM_DIR"/S18_marked.bam "$BAM_DIR"/S325_marked.bam "$BAM_DIR"/S344_marked.bam "$BAM_DIR"/S39_marked.bam "$BAM_DIR"/S19_marked.bam "$BAM_DIR"/S326_marked.bam "$BAM_DIR"/S345_marked.bam "$BAM_DIR"/S3_marked.bam "$BAM_DIR"/S1_marked.bam "$BAM_DIR"/S327_marked.bam "$BAM_DIR"/S346_marked.bam "$BAM_DIR"/S40_marked.bam "$BAM_DIR"/S20_marked.bam "$BAM_DIR"/S328_marked.bam "$BAM_DIR"/S347_marked.bam "$BAM_DIR"/S4_marked.bam "$BAM_DIR"/S21_marked.bam "$BAM_DIR"/S329_marked.bam "$BAM_DIR"/S348_marked.bam "$BAM_DIR"/S5_marked.bam "$BAM_DIR"/S22_marked.bam "$BAM_DIR"/S32_marked.bam "$BAM_DIR"/S349_marked.bam "$BAM_DIR"/S6_marked.bam "$BAM_DIR"/S23_marked.bam "$BAM_DIR"/S330_marked.bam "$BAM_DIR"/S34_marked.bam "$BAM_DIR"/S7_marked.bam "$BAM_DIR"/S24_marked.bam "$BAM_DIR"/S331_marked.bam "$BAM_DIR"/S350_marked.bam "$BAM_DIR"/S8_marked.bam "$BAM_DIR"/S25_marked.bam "$BAM_DIR"/S332_marked.bam "$BAM_DIR"/S351_marked.bam "$BAM_DIR"/S9_marked.bam "$BAM_DIR"/S26_marked.bam "$BAM_DIR"/S333_marked.bam "$BAM_DIR"/S352_marked.bam "$BAM_DIR"/S27_marked.bam "$BAM_DIR"/S334_marked.bam "$BAM_DIR"/S353_marked.bam "$BAM_DIR"/S28_marked.bam "$BAM_DIR"/S335_marked.bam "$BAM_DIR"/S354_marked.bam "$BAM_DIR"/S29_marked.bam "$BAM_DIR"/S336_marked.bam "$BAM_DIR"/S355_marked.bam
done
```

Compressed and indexed WhatsHap output.
```bash
module load bcftools
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)
for NAME in "${VCFLIST[@]}"; do bgzip "$NAME"_hapmarked.vcf; bcftools index "$NAME"_hapmarked.vcf.gz; done
```

Reformatted genetic map into HapMap format and then ran SHAPEIT4.
```bash
# Ran in UFRC queue system; see shapeit.job for more details.
# Resources Used: 350 Mb, 5 min

module load shapeit4/4.2.2
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
declare -a VCFLIST=(01 02 03 04 05 06 07 08 09 10 11)

for NUM in "${VCFLIST[@]}"
do
    shapeit4 --input chr"$NUM"_hapmarked.vcf.gz --map "$MAP_DIR"/1060_LH_F2_manual_copy_Chr"$NUM".gmap --region Chr"$NUM" --output chr"$NUM"_phased.vcf --thread 12
done
```

## E. globulus LD
### Whole-genome estimates
### Minor Allele Count = 1
Used [`emeraLD`](https://github.com/statgen/emeraLD/tree/master) to calculate LD in 100kb non-overlapping sliding windows across all chromosomes.

Used `tabix` to index phased variants
```bash
module load htslib/1.15
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
declare -a VCFLIST=(01 02 03 04 05 06 07 08 09 10 11)

cd "$WDIR"
for NUM in "${VCFLIST[@]}"
do
    # bgzip chr"$NUM"_hapmarked.vcf
    bgzip chr"$NUM"_phased.vcf
    tabix chr"$NUM"_phased.vcf.gz
done
```

Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.
```bash
# Ran on UFRC queue system; see glob_genwide_ld_mac01.job for more details.
# Resources used: 70 Mb, 2 hrs

module load htslib/1.15
module load emerald/0.1

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/genome_wide/mac01"

emeraLD --in "$VCF_DIR"/chr01_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr01_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr02_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr02_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr03_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr03_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr04_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr04_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr05_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr05_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr06_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr06_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr07_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr07_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr08_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr08_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr09_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr09_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr10_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr10_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr11_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$WDIR"/glob_Chr11_ld.txt.gz
```
Averaged r2 values for each distance between variants using a custom `python` script.
```bash
# Done in UFRC queue system; see glob_genwide_mac01_avg.job
# Resources used: 21 Gb, 4 hrs

module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/genome_wide/mac01"
declare -a NAMELIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPT_DIR"/average_r2.py -i "$TAB_DIR"/glob_"$NAME"_ld.txt.gz -o glob_"$NAME"_ld_avg.csv
done
```

Fitted LD curve to pairwise r2 values versus marker distance in `R` using equation from Hill and Weir 1988 and custom `R` script.

```bash
# Done in UFRC queue system; see genwide_mac01_fit.job for more details.
# Resources used: 630 Mb, 7 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/genome_wide"

Rscript "$SCRIPT_DIR"/fit_rsq_curve.r -f chr_r2_tab_list.txt -o genwide_ld.txt -m 100 -r 0.2 -g TRUE
```

| Chromosome | LD r2 = 0.2 (bp) | C        |
| ---------- | ---------------- | -------- |
| Chr01      | 1983             | 0.001761 |
| Chr02      | 1912             | 0.001826 |
| Chr03      | 1552             | 0.002249 |
| Chr04      | 1572             | 0.002221 |
| Chr05      | 1207             | 0.002893 |
| Chr06      | 1814             | 0.001924 |
| Chr07      | 902              | 0.003870 |
| Chr08      | 1832             | 0.001906 |
| Chr09      | 2504             | 0.001394 |
| Chr10      | 2689             | 0.001299 |
| Chr11      | 1332             | 0.002621 |
| AVERAGE    | 1754.45          | 0.001977 |

### Minor Allele Count = 2
Equivalent to minor allele frequency = 0.05. Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.

```bash
# Done in UFRC queue system; see glob_genwide_mac02_ld.job for more details.
# Resources used: 65 Mb, 1.5 hrs

module purge
module load htslib/1.15
module load emerald/0.1

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"

emeraLD --in "$VCF_DIR"/chr01_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr01_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr02_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr02_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr03_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr03_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr04_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr04_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr05_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr05_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr06_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr06_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr07_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr07_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr08_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr08_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr09_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr09_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr10_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr10_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr11_phased.vcf.gz --phase --mac 2 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac02/glob_Chr11_ld.txt.gz
```

Averaged r2 values for each distance between variants using a custom `python` script.
```bash
# Done in UFRC queue system; see glob_genwide_mac02_avg.job for more details.
# Resources used: 17 Gb, 3 hrs

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/genome_wide/mac02"
declare -a NAMELIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPT_DIR"/average_r2.py -i "$TAB_DIR"/glob_"$NAME"_ld.txt.gz -o glob_"$NAME"_ld_avg.csv
done
```

Fitted LD curve to pairwise r2 values versus marker distance in `R` using equation from Hill and Weir 1988 and custom `R` script.

```bash
# Done in UFRC queue system; see glob_genwide_mac02_fit.job for more details.
# Resources used: 700 Mb, 5 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/genome_wide"

Rscript "$SCRIPT_DIR"/fit_rsq_curve.r -f chr_r2_tab_list.txt -o genwide_ld.txt -m 100 -r 0.2 -g TRUE
```

| Chromosome | LD r2 = 0.2 (bp) | C        |
| ---------- | ---------------- | -------- |
| Chr01      | 2261             | 0.001544 |
| Chr02      | 2450             | 0.001425 |
| Chr03      | 1964             | 0.001778 |
| Chr04      | 1964             | 0.001778 |
| Chr05      | 1449             | 0.002409 |
| Chr06      | 2220             | 0.001573 |
| Chr07      | 1125             | 0.003104 |
| Chr08      | 2224             | 0.001570 |
| Chr09      | 3120             | 0.001119 |
| Chr10      | 3331             | 0.001048 |
| Chr11      | 1674             | 0.002086 |
| AVERAGE    | 2162             | 0.001767 |

### Sliding window estimates
Generated jobs to run emeraLD across sliding windows, using MAC = 1 filtering and a window size of 100,000 bp (default values for job-making script).

```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/windows"

python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr01_phased.vcf.gz -c Chr01 -s 42219553 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr01/glob_chr01_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr02_phased.vcf.gz -c Chr02 -s 50828380 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr02/glob_chr02_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr03_phased.vcf.gz -c Chr03 -s 65547241 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr03/glob_chr03_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr04_phased.vcf.gz -c Chr04 -s 38599333 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr04/glob_chr04_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr05_phased.vcf.gz -c Chr05 -s 62919971 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr05/glob_chr05_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr06_phased.vcf.gz -c Chr06 -s 52140660 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr06/glob_chr06_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr07_phased.vcf.gz -c Chr07 -s 54252628 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr07/glob_chr07_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr08_phased.vcf.gz -c Chr08 -s 70214608 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr08/glob_chr08_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr09_phased.vcf.gz -c Chr09 -s 38300324 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr09/glob_chr09_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr10_phased.vcf.gz -c Chr10 -s 38722660 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr10/glob_chr10_wins_ld.job -t glob
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr11_phased.vcf.gz -c Chr11 -s 42056460 -i "$LIST_DIR"/Eglobulus_MR.txt -o "$WDIR"/chr11/glob_chr11_wins_ld.job -t glob
```

Not displaying example of jobfiles here since they are quite large; see chrXX_ld_mac01.job in the jobfiles directory. Ran emeraLD on each window separately and then averaged r2 for each distance class within individual windows. Each job used about 15 Mb RAM and took ~2 minutes.

```bash
# Done in UFRC queue system; see glob_window_mac01_avg.job for more details.
# Resources used: 130 Mb, 2 hrs

module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/windows"
BASE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/windows"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    cd "$BASE_DIR"/"$NAME"
    ls *.txt.gz | while read FILE
    do
        python "$SCRIPT_DIR"/average_r2.py -i "$BASE_DIR"/"$NAME"/"$FILE" -o "$WDIR"/"$NAME"/"$FILE"_avg.csv
    done
done
```

Created file lists for curve fitting.
```bash
declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)
for NAME in "${NAMELIST[@]}"; do ls /blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/windows/"$NAME"/*.csv > glob_"$NAME"_win_ld_files.txt; done
```

Fit LD decay curves for each window and interpolated distance at which r^2 = 0.2 using a custom `R` script.

```bash
# Done in UFRC queue system; see glob_window_mac01_cfit.job for more details.
# Resources used: 320 Mb, 20 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/curve_fit"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    Rscript "$SCRIPT_DIR"/fit_rsq_curve.r -f "$NAME"_win_ld_files.txt -o glob_"$NAME"_windows_ld.txt -m 40 -r 0.2 -g FALSE
done
```

Plotted and curve fitted LD across each chromosome in `R`.
```R
loess_graph <- function(chr_ld_tab, smoothing, plot_title){
    chr_ld_tab <- chr_ld_tab[!is.na(chr_ld_tab$ld),]
    temp <- unlist(lapply(strsplit(chr_ld_tab$file_name, "_"), function(x) x[2]))
    chr_ld_tab$pos <- as.numeric(unlist(lapply(strsplit(temp, "-"), function(x) x[1])))
    chr_ld_tab <- chr_ld_tab[order(chr_ld_tab$pos),]
    chr_loess <- loess(ld ~ pos, data=chr_ld_tab, span=smoothing)
    chr_pred <- predict(chr_loess)
    plot(chr_ld_tab$pos, chr_ld_tab$ld, main = plot_title, xlab = "bp", ylab = "ld", pch = 19)
    lines(chr_pred, x = chr_loess$x, col = "blue")
    chr_ld_tab_sorted <- chr_ld_tab[order(chr_ld_tab$ld), ]
    outlier_tab <- chr_ld_tab_sorted[c(round(nrow(chr_ld_tab_sorted) * 0.95):nrow(chr_ld_tab_sorted)),]
    outlier_tab <- outlier_tab[order(outlier_tab$pos), ]
    return(outlier_tab)
}

# Chr01
chr01_ld_tab <- read.table("glob_chr01_windows_ld.txt", header = TRUE, sep = "\t")
chr01_outliers <- loess_graph(chr01_ld_tab, 0.1, "E.glob Chr01 LD MAC = 0.1")
# Chr02
chr02_ld_tab <- read.table("glob_chr02_windows_ld.txt", header = TRUE, sep = "\t")
chr02_outliers <- loess_graph(chr02_ld_tab, 0.1, "E.glob Chr02 LD MAC = 0.1")
# Chr03
chr03_ld_tab <- read.table("glob_chr03_windows_ld.txt", header = TRUE, sep = "\t")
chr03_outliers <- loess_graph(chr03_ld_tab, 0.1, "E.glob Chr03 LD MAC = 0.1")
# Chr04
chr04_ld_tab <- read.table("glob_chr04_windows_ld.txt", header = TRUE, sep = "\t")
chr04_outliers <- loess_graph(chr04_ld_tab, 0.1, "E.glob Chr04 LD MAC = 0.1")
# Chr05
chr05_ld_tab <- read.table("glob_chr05_windows_ld.txt", header = TRUE, sep = "\t")
chr05_outliers <- loess_graph(chr05_ld_tab, 0.1, "E.glob Chr05 LD MAC = 0.1")
# Chr06
chr06_ld_tab <- read.table("glob_chr06_windows_ld.txt", header = TRUE, sep = "\t")
chr06_outliers <- loess_graph(chr06_ld_tab, 0.1, "E.glob Chr06 LD MAC = 0.1")
# Chr07
chr07_ld_tab <- read.table("glob_chr07_windows_ld.txt", header = TRUE, sep = "\t")
chr07_outliers <- loess_graph(chr07_ld_tab, 0.1, "E.glob Chr07 LD MAC = 0.1")
# Chr08
chr08_ld_tab <- read.table("glob_chr08_windows_ld.txt", header = TRUE, sep = "\t")
chr08_outliers <- loess_graph(chr08_ld_tab, 0.1, "E.glob Chr08 LD MAC = 0.1")
# Chr09
chr09_ld_tab <- read.table("glob_chr09_windows_ld.txt", header = TRUE, sep = "\t")
chr09_outliers <- loess_graph(chr09_ld_tab, 0.1, "E.glob Chr09 LD MAC = 0.1")
# Chr10
chr10_ld_tab <- read.table("glob_chr10_windows_ld.txt", header = TRUE, sep = "\t")
chr10_outliers <- loess_graph(chr10_ld_tab, 0.1, "E.glob Chr10 LD MAC = 0.1")
# Chr11
chr11_ld_tab <- read.table("glob_chr11_windows_ld.txt", header = TRUE, sep = "\t")
chr11_outliers <- loess_graph(chr11_ld_tab, 0.1, "E.glob Chr11 LD MAC = 0.1")
```

Intervals of interest (LD outliers):
| Chr   | Interval            | LD         |
| ----- | ------------------- | ---------- |
| Chr01 | 29500000 - 29699999 | 11430      |
| Chr02 | 39600000 - 40099999 | 7841       |
| Chr02 | 46700000 - 46999999 | 12732      |
| Chr03 | 29800000 - 29999999 | 15773      |
| Chr04 |  9700000 -  9899999 | 10338      |
| Chr04 | 15100000 - 15899999 | 11680      |
| Chr05 | 30100000 - 30399999 | 7708       |
| Chr05 | 61800000 - 61999999 | 18512      |
| Chr06 | 43500000 - 43699999 | 8773       |
| Chr07 | 22500000 - 22799999 | 8116       |
| Chr07 | 50400000 - 50699999 | 18108      |
| Chr08 | 11200000 - 11399999 | 12763      |
| Chr08 | 36000000 - 36699999 | 20992      |
| Chr09 | 22600000 - 22799999 | 42340      |
| Chr09 | 24300000 - 24999999 | 19976      |
| Chr09 | 27700000 - 27999999 | 24911      |
| Chr10 |  4000000 -  4299999 | 43987      |
| Chr10 | 23800000 - 23999999 | 19112      |
| Chr11 | 24900000 - 25499999 | 11549      |

Consolidated window output files per chromosome into one genome-wide.
```bash
declare -a CHRLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)
head -n 1 glob_chr01_windows_ld.txt > glob_all_windows_ld.txt

for NAME in "${CHRLIST[@]}"
do
    tail -n +2 glob_"$NAME"_windows_ld.txt >> glob_all_windows_ld.txt
done
```

## E. cordata LD
### Genome-wide MAC = 1
Equivalent to minor allele frequency = 0.05. Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.

```bash
# Done in UFRC queue system; see cord_genwide_mac01_ld.job for more details.
# Resources used: 60 Mb, 1.5 hrs

module purge
module load htslib/1.15
module load emerald/0.1

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/emerald/genome_wide"

emeraLD --in "$VCF_DIR"/chr01_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr01_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr02_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr02_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr03_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr03_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr04_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr04_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr05_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr05_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr06_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr06_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr07_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr07_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr08_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr08_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr09_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr09_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr10_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr10_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr11_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$WDIR"/cord_Chr11_ld.txt.gz
```

Averaged r2 values for each distance between variants using a custom `python` script.
```bash
# Done in UFRC queue system; see cord_genwide_mac01_avg.job for more details.
# Resources used: 13 Gb, 2 hrs

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/emerald/genome_wide"
declare -a NAMELIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPT_DIR"/average_r2.py "$TAB_DIR"/cord_"$NAME"_ld.txt.gz cord_"$NAME"_ld_avg.csv
done
```

Fitted LD curve to pairwise r2 values versus marker distance in `R` using equation from Hill and Weir 1988 and custom `R` script.

```bash
# Done in UFRC queue system; see cord_genwide_mac01_fit.job for more details.
# Resources used: 620 Mb, 5 min

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/emerald/genome_wide"
declare -a NAMELIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPT_DIR"/average_r2.py "$TAB_DIR"/cord_"$NAME"_ld.txt.gz cord_"$NAME"_ld_avg.csv
done
```

| Chromosome | LD r2 = 0.2 (bp) | C         |
| ---------- | ---------------- | --------- |
| Chr01      | 20501            | 0.0001703 |
| Chr02      | 18654            | 0.0001871 |
| Chr03      | 17819            | 0.0001959 |
| Chr04      | 15812            | 0.0002208 |
| Chr05      | 16466            | 0.0002120 |
| Chr06      | 16314            | 0.0002140 |
| Chr07      | 19961            | 0.0001749 |
| Chr08      | 14610            | 0.0002390 |
| Chr09      | 17393            | 0.0002007 |
| Chr10      | 15322            | 0.0002278 |
| Chr11      | 16905            | 0.0002065 |
| AVERAGE    | 17,250.63        | 0.0002203 |

### Sliding window estimates
Generated jobs to run emeraLD across sliding windows, using MAC = 1 filtering and sliding window size = 100kb (defaults for script).
```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/emerald/windows"

python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr01_phased.vcf.gz -c Chr01 -s 42219553 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr01/cord_chr01_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr02_phased.vcf.gz -c Chr02 -s 50828380 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr02/cord_chr02_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr03_phased.vcf.gz -c Chr03 -s 65547241 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr03/cord_chr03_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr04_phased.vcf.gz -c Chr04 -s 38599333 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr04/cord_chr04_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr05_phased.vcf.gz -c Chr05 -s 62919971 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr05/cord_chr05_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr06_phased.vcf.gz -c Chr06 -s 52140660 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr06/cord_chr06_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr07_phased.vcf.gz -c Chr07 -s 54252628 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr07/cord_chr07_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr08_phased.vcf.gz -c Chr08 -s 70214608 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr08/cord_chr08_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr09_phased.vcf.gz -c Chr09 -s 38300324 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr09/cord_chr09_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr10_phased.vcf.gz -c Chr10 -s 38722660 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr10/cord_chr10_wins_ld.job -t cord
python "$SCRIPT_DIR"/make_emerald_job.py -v "$VCF_DIR"/chr11_phased.vcf.gz -c Chr11 -s 42056460 -i "$LIST_DIR"/Ecordata.txt -o "$WDIR"/chr11/cord_chr11_wins_ld.job -t cord
```
Not displaying example of jobfiles here since they are quite large; see chrXX_ld_mac01.job in the jobfiles directory. Ran emeraLD on each window separately and then averaged r2 for each distance class within individual windows. All jobs ran with a maximum of 15 Mb RAM and 2 minutes.

```bash
# Done in UFRC queue system; see cord_window_mac01_avg.job for more details.
# Resources used: 95 Mb, 2 hrs

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/avg_r2"
BASE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/emerald"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    cd "$BASE_DIR"/"$NAME"
    ls *.txt.gz | while read FILE
    do
        python "$SCRIPT_DIR"/average_r2.py -i "$BASE_DIR"/"$NAME"/"$FILE" -o "$WDIR"/"$NAME"/"$FILE"_avg.csv
    done
done
```

Created file lists for curve fitting.
```bash
declare -a CHRLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)
for CHR in "${CHRLIST[@]}"; do ls /blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/avg_r2/windows/"$CHR"/*.csv > cord_"$CHR"_win_ld_files.txt; done
```

Fit LD decay curves for each window and interpolated distance at which r^2 = 0.2 using a custom `R` script.

```bash
# Done in UFRC queue system; see glob_window_mac01_cfit.job for more details.
# Resources used: 200 Mb, 15 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord/curve_fit"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    Rscript "$SCRIPT_DIR"/fit_rsq_curve.r -f cord_"$NAME"_win_ld_files.txt -o cord_"$NAME"_windows_ld.txt -m 40 -r 0.2 -g FALSE
done
```

Plotted and curve fitted LD across each chromosome in `R`.
```R
loess_graph <- function(chr_ld_tab, smoothing, plot_title){
    chr_ld_tab <- chr_ld_tab[!is.na(chr_ld_tab$ld),]
    temp <- unlist(lapply(strsplit(chr_ld_tab$file_name, "_"), function(x) x[2]))
    chr_ld_tab$pos <- as.numeric(unlist(lapply(strsplit(temp, "-"), function(x) x[1])))
    chr_ld_tab <- chr_ld_tab[order(chr_ld_tab$pos),]
    chr_loess <- loess(ld ~ pos, data=chr_ld_tab, span=smoothing)
    chr_pred <- predict(chr_loess)
    plot(chr_ld_tab$pos, chr_ld_tab$ld, main = plot_title, xlab = "bp", ylab = "ld", pch = 19)
    lines(chr_pred, x = chr_loess$x, col = "blue")
    chr_ld_tab_sorted <- chr_ld_tab[order(chr_ld_tab$ld), ]
    outlier_tab <- chr_ld_tab_sorted[c(round(nrow(chr_ld_tab_sorted) * 0.95):nrow(chr_ld_tab_sorted)),]
    outlier_tab <- outlier_tab[order(outlier_tab$pos), ]
    return(outlier_tab)
}

# Chr01
chr01_ld_tab <- read.table("cord_chr01_windows_ld.txt", header = TRUE, sep = "\t")
chr01_outliers <- loess_graph(chr01_ld_tab, 0.1, "E.cord Chr01 LD MAC = 0.1")
# Chr02
chr02_ld_tab <- read.table("cord_chr02_windows_ld.txt", header = TRUE, sep = "\t")
chr02_outliers <- loess_graph(chr02_ld_tab, 0.1, "E.cord Chr02 LD MAC = 0.1")
# Chr03
chr03_ld_tab <- read.table("cord_chr03_windows_ld.txt", header = TRUE, sep = "\t")
chr03_outliers <- loess_graph(chr03_ld_tab, 0.1, "E.cord Chr03 LD MAC = 0.1")
# Chr04
chr04_ld_tab <- read.table("cord_chr04_windows_ld.txt", header = TRUE, sep = "\t")
chr04_outliers <- loess_graph(chr04_ld_tab, 0.1, "E.cord Chr04 LD MAC = 0.1")
# Chr05
chr05_ld_tab <- read.table("cord_chr05_windows_ld.txt", header = TRUE, sep = "\t")
chr05_outliers <- loess_graph(chr05_ld_tab, 0.1, "E.cord Chr05 LD MAC = 0.1")
# Chr06
chr06_ld_tab <- read.table("cord_chr06_windows_ld.txt", header = TRUE, sep = "\t")
chr06_outliers <- loess_graph(chr06_ld_tab, 0.1, "E.cord Chr06 LD MAC = 0.1")
# Chr07
chr07_ld_tab <- read.table("cord_chr07_windows_ld.txt", header = TRUE, sep = "\t")
chr07_outliers <- loess_graph(chr07_ld_tab, 0.1, "E.cord Chr07 LD MAC = 0.1")
# Chr08
chr08_ld_tab <- read.table("cord_chr08_windows_ld.txt", header = TRUE, sep = "\t")
chr08_outliers <- loess_graph(chr08_ld_tab, 0.1, "E.cord Chr08 LD MAC = 0.1")
# Chr09
chr09_ld_tab <- read.table("cord_chr09_windows_ld.txt", header = TRUE, sep = "\t")
chr09_outliers <- loess_graph(chr09_ld_tab, 0.1, "E.cord Chr09 LD MAC = 0.1")
# Chr10
chr10_ld_tab <- read.table("cord_chr10_windows_ld.txt", header = TRUE, sep = "\t")
chr10_outliers <- loess_graph(chr10_ld_tab, 0.1, "E.cord Chr10 LD MAC = 0.1")
# Chr11
chr11_ld_tab <- read.table("cord_chr11_windows_ld.txt", header = TRUE, sep = "\t")
chr11_outliers <- loess_graph(chr11_ld_tab, 0.1, "E.cord Chr11 LD MAC = 0.1")
```

Consolidated window output files per chromosome into one genome-wide.
```bash
declare -a CHRLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)
head -n 1 cord_chr01_windows_ld.txt > cord_all_windows_ld.txt

for NAME in "${CHRLIST[@]}"
do
    tail -n +2 cord_"$NAME"_windows_ld.txt >> cord_all_windows_ld.txt
done
```