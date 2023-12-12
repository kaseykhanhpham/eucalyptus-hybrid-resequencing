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
# Resources Used: ??

module load whatshap/1.1
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
BAM_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    whatshap phase -o "$NAME"_hapmarked.vcf --reference "$REF_DIR"/EGLOB-X46.v1.0.fa --tag PS "$VCF_DIR"/"$NAME"_fil.vcf.gz "$BAM_DIR"/S10_marked.bam "$BAM_DIR"/S2_marked.bam "$BAM_DIR"/S337_marked.bam "$BAM_DIR"/S356_marked.bam "$BAM_DIR"/S11_marked.bam "$BAM_DIR"/S30_marked.bam "$BAM_DIR"/S338_marked.bam "$BAM_DIR"/S357_marked.bam "$BAM_DIR"/S12_marked.bam "$BAM_DIR"/S31_marked.bam "$BAM_DIR"/S339_marked.bam "$BAM_DIR"/S358_marked.bam "$BAM_DIR"/S13_marked.bam "$BAM_DIR"/S320_marked.bam "$BAM_DIR"/S33_marked.bam "$BAM_DIR"/S359_marked.bam "$BAM_DIR"/S14_marked.bam "$BAM_DIR"/S321_marked.bam "$BAM_DIR"/S340_marked.bam "$BAM_DIR"/S35_marked.bam "$BAM_DIR"/S15_marked.bam "$BAM_DIR"/S322_marked.bam "$BAM_DIR"/S341_marked.bam "$BAM_DIR"/S36_marked.bam "$BAM_DIR"/S16_marked.bam "$BAM_DIR"/S323_marked.bam "$BAM_DIR"/S342_marked.bam "$BAM_DIR"/S37_marked.bam "$BAM_DIR"/S17_marked.bam "$BAM_DIR"/S324_marked.bam "$BAM_DIR"/S343_marked.bam "$BAM_DIR"/S38_marked.bam "$BAM_DIR"/S18_marked.bam "$BAM_DIR"/S325_marked.bam "$BAM_DIR"/S344_marked.bam "$BAM_DIR"/S39_marked.bam "$BAM_DIR"/S19_marked.bam "$BAM_DIR"/S326_marked.bam "$BAM_DIR"/S345_marked.bam "$BAM_DIR"/S3_marked.bam "$BAM_DIR"/S1_marked.bam "$BAM_DIR"/S327_marked.bam "$BAM_DIR"/S346_marked.bam "$BAM_DIR"/S40_marked.bam "$BAM_DIR"/S20_marked.bam "$BAM_DIR"/S328_marked.bam "$BAM_DIR"/S347_marked.bam "$BAM_DIR"/S4_marked.bam "$BAM_DIR"/S21_marked.bam "$BAM_DIR"/S329_marked.bam "$BAM_DIR"/S348_marked.bam "$BAM_DIR"/S5_marked.bam "$BAM_DIR"/S22_marked.bam "$BAM_DIR"/S32_marked.bam "$BAM_DIR"/S349_marked.bam "$BAM_DIR"/S6_marked.bam "$BAM_DIR"/S23_marked.bam "$BAM_DIR"/S330_marked.bam "$BAM_DIR"/S34_marked.bam "$BAM_DIR"/S7_marked.bam "$BAM_DIR"/S24_marked.bam "$BAM_DIR"/S331_marked.bam "$BAM_DIR"/S350_marked.bam "$BAM_DIR"/S8_marked.bam "$BAM_DIR"/S25_marked.bam "$BAM_DIR"/S332_marked.bam "$BAM_DIR"/S351_marked.bam "$BAM_DIR"/S9_marked.bam "$BAM_DIR"/S26_marked.bam "$BAM_DIR"/S333_marked.bam "$BAM_DIR"/S352_marked.bam "$BAM_DIR"/SRR10339635_marked.bam "$BAM_DIR"/S27_marked.bam "$BAM_DIR"/S334_marked.bam "$BAM_DIR"/S353_marked.bam "$BAM_DIR"/S28_marked.bam "$BAM_DIR"/S335_marked.bam "$BAM_DIR"/S354_marked.bam "$BAM_DIR"/S29_marked.bam "$BAM_DIR"/S336_marked.bam "$BAM_DIR"/S355_marked.bam
done
```

Compressed and indexed WhatsHap output.
```bash
module load bcftools
for NAME in "${VCFLIST[@]}"; do bgzip "$NAME"_hapmarked.vcf; bcftools index "$NAME"_hapmarked.vcf.gz; done
```

Reformatted genetic map into HapMap format and then ran SHAPEIT4.
```bash
# Ran in UFRC queue system; see shapeit.job for more details.
# Resources Used: 330 Mb, 4 min

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
    bgzip chr"$NUM"_hapmarked.vcf
    bgzip chr"$NUM"_phased.vcf
    tabix chr"$NUM"_phased.vcf.gz
done
```

Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.
```bash
# Ran on UFRC queue system; see glob_genwide_ld_mac01.job for more details.
# Resources used: 50 Mb, 1 hr

module purge
module load htslib/1.15
module load emerald/0.1

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"

emeraLD --in "$VCF_DIR"/chr01_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr01_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr02_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr02_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr03_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr03_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr04_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr04_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr05_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr05_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr06_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr06_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr07_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr07_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr08_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr08_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr09_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr09_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr10_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr10_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr11_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Eglobulus_MR.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/glob_Chr11_ld.txt.gz
```
Averaged r2 values for each distance between variants using a custom `python` script.
```

```

Fitted LD curve to pairwise r2 values versus marker distance in `R` using equation from Hill and Weir 1988 and custom `R` script.

```bash
# Done in UFRC queue system; see genwide_mac01_fit.job for more details.
# Resources used: 500 Mb, 5 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/genome_wide"

Rscript "$SCRIPT_DIR"/fit_rsq_curve.r chr_r2_tab_list.txt genwide_ld.txt 100 0.2 TRUE
```

| Chromosome | LD r2 = 0.2 (bp) | C         |
| ---------- | ---------------- | --------- |
| Chr01      | 1926             | 0.001813  |
| Chr02      | 1795             | 0.001945  |
| Chr03      | 1363             | 0.002562  |
| Chr04      | 1474             | 0.002369  |
| Chr05      | 1184             | 0.002948  |
| Chr06      | 1647             | 0.002119  |
| Chr07      | 932              | 0.003747  |
| Chr08      | 1659             | 0.002104  |
| Chr09      | 2278             | 0.001533  |
| Chr10      | 2412             | 0.001447  |
| Chr11      | 1251             | 0.002792  |
| AVERAGE    | 1629             | 0.002307  |

### Minor Allele Count = 2
Equivalent to minor allele frequency = 0.05. Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.

```bash
# Done in UFRC queue system; see glob_genwide_mac02_ld.job for more details.
# Resources used: 45 Mb, 30 min

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
# Resources used: 7 Gb, 2 hrs

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald/genome_wide/mac02"
declare -a NAMELIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPT_DIR"/average_r2.py "$TAB_DIR"/glob_"$NAME"_ld.txt.gz glob_"$NAME"_ld_avg.csv
done
```

Fitted LD curve to pairwise r2 values versus marker distance in `R` using equation from Hill and Weir 1988 and custom `R` script.

```bash
# Done in UFRC queue system; see glob_genwide_mac02_fit.job for more details.
# Resources used: 700 Mb, 7 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
TAB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2/genome_wide"

Rscript "$SCRIPT_DIR"/fit_rsq_curve.r chr_r2_tab_list.txt genwide_ld.txt 100 0.2 TRUE
```

| Chromosome | LD r2 = 0.2 (bp) | C         |
| ---------- | ---------------- | --------- |
| Chr01      | 2228             | 0.001567  |
| Chr02      | 2331             | 0.001498  |
| Chr03      | 1748             | 0.001997  |
| Chr04      | 1864             | 0.001873  |
| Chr05      | 1430             | 0.002441  |
| Chr06      | 1982             | 0.001761  |
| Chr07      | 1173             | 0.002978  |
| Chr08      | 2024             | 0.001725  |
| Chr09      | 2844             | 0.001228  |
| Chr10      | 3011             | 0.001159  |
| Chr11      | 1582             | 0.002207  |
| AVERAGE    | 2020             | 0.001858  |

### Sliding window estimates
Generated jobs to run emeraLD across windows in each chromosome
```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr01_phased.vcf.gz Chr01 42219553 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr01/mac01/chr01_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr02_phased.vcf.gz Chr02 50828380 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr02/mac01/chr02_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr03_phased.vcf.gz Chr03 65547241 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr03/mac01/chr03_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr04_phased.vcf.gz Chr04 38599333 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr04/mac01/chr04_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr05_phased.vcf.gz Chr05 62919971 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr05/mac01/chr05_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr06_phased.vcf.gz Chr06 52140660 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr06/mac01/chr06_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr07_phased.vcf.gz Chr07 54252628 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr07/mac01/chr07_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr08_phased.vcf.gz Chr08 70214608 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr08/mac01/chr08_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr09_phased.vcf.gz Chr09 38300324 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr09/mac01/chr09_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr10_phased.vcf.gz Chr10 38722660 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr10/mac01/chr10_ld_mac01.job
python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr11_phased.vcf.gz Chr11 42056460 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr11/mac01/chr11_ld_mac01.job
```

Not displaying example of jobfiles here since they are quite large; see chrXX_ld_mac01.job in the jobfiles directory. Ran emeraLD on each window separately and then averaged r2 for each distance class within individual windows.

```bash
# Done in UFRC queue system; see glob_window_mac01_avg.job for more details.
# Resources used: 

module purge
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/avg_r2"
BASE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/emerald"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    cd "$BASE_DIR"/"$NAME"/mac01
    ls *.txt.gz | while read FILE
    do
        python "$SCRIPT_DIR"/average_r2.py "$BASE_DIR"/"$NAME"/mac01/"$FILE" "$WDIR"/"$NAME"/mac01/"$FILE"_avg.csv
    done
done
```

Fit LD decay curves for each window and interpolated distance at which r^2 = 0.2 using a custom `R` script.

```bash
# Done in UFRC queue system; see glob_window_mac01_cfit.job for more details.
# Resources used: 200 Mb, 15 min

module purge
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/curve_fit"

declare -a NAMELIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11)

for NAME in "${NAMELIST[@]}"
do
    Rscript "$SCRIPT_DIR"/fit_rsq_curve.r "$NAME"_win_ld_files.txt glob_"$NAME"_windows_ld.txt 40 0.2 FALSE
done
```

Plotted and curve fitted LD across each chromosome in `R`.
```R
# Chr01
chr01_ld_tab <- read.table("glob_chr01_windows_ld.txt", header = TRUE, sep = "\t")
temp <- unlist(lapply(strsplit(chr01_ld_tab$file_name, "_"), function(x) x[2]))
chr01_ld_tab$pos <- unlist(lapply(strsplit(temp, "-"), function(x) x[1]))
chr01_loess <- loess(ld ~ pos, data=chr01_ld_tab, span=0.1)
chr01_pred <- predict(chr01_loess)
plot(chr01_ld_tab$pos, chr01_ld_tab$ld, main = "Chr01 LD MAC=1", xlab = "bp", ylab = "ld")
lines(chr01_pred, x = chr01_ld_tab$pos, col = "blue")
```

## E. cordata LD
### Genome-wide MAC = 1
Equivalent to minor allele frequency = 0.05. Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.

```bash
# Done in UFRC queue system; see cord_genwide_mac01_ld.job for more details.
# Resources used: 41 Mb, 30 min

module purge
module load htslib/1.15
module load emerald/0.1

VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/cord"

emeraLD --in "$VCF_DIR"/chr01_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr01_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr02_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr02_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr03_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr03_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr04_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr04_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr05_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr05_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr06_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr06_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr07_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr07_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr08_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr08_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr09_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr09_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr10_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/mac01/cord_Chr10_ld.txt.gz
emeraLD --in "$VCF_DIR"/chr11_phased.vcf.gz --phase --mac 1 --include "$LIST_DIR"/Ecordata.txt --stdout | bgzip -c > "$LIST_DIR"/emerald/genome_wide/cord_Chr11_ld.txt.gz
```

Averaged r2 values for each distance between variants using a custom `python` script.
```bash
# Done in UFRC queue system; see cord_genwide_mac01_avg.job for more details.
# Resources used: 6 Gb, 1 hr

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
# Resources used: 600 Mb, 6 min

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
| Chr01      | 19457            | 0.0001794 |
| Chr02      | 19248            | 0.0001814 |
| Chr03      | 18195            | 0.0001919 |
| Chr04      | 15220            | 0.0002294 |
| Chr05      | 17203            | 0.0002030 |
| Chr06      | 15763            | 0.0002215 |
| Chr07      | 19716            | 0.0001771 |
| Chr08      | 15657            | 0.0002230 |
| Chr09      | 17680            | 0.0001975 |
| Chr10      | 15056            | 0.0002319 |
| Chr11      | 16747            | 0.0002085 |
| AVERAGE    | 17267            | 0.0002040 |