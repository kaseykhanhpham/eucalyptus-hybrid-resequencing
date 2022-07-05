# Linkage Disequilibrium Calculation

[Meier and Ravinet's tutorial](https://speciationgenomics.github.io/ld_decay/) referenced extensively.

## MAF = 0.025

### Long Range

**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.025_long.job for more information.
# Resources: 10 Gb, 30 min

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld"

while read CHR
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.025 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 100000 --ld-window-kb 1000 -ld-window-r2 0 --make-bed  --vcf-half-call m --thin 0.5 --out "$CHR" --threads 12
done < "$LIST_DIR"/chromosome_list.txt
```

I realized that it wouldn't be worth it to use the unanchored contigs for ld calculation because they contain too few SNPs, so I deleted those files.

**Process plink outputs:**

Wrote python script to average pairwise R2 values for each distance between SNPs and ran for each chromosome.
```bash
# Run via UFRC queue, see avg_r2_maf0.025_long.job for more details.
# Resources: 45 Gb, 10 hrs

module load python/3.8

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/maf0.025/long_range"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for NAME in "${CHROMOSOMES[@]}"
do
        python "$SCRIPTS_DIR"/average_r2.py  "$WDIR"/"$NAME".ld.gz "$WDIR"/"$NAME"_r2.csv
done
```

**Plot plink output:**

Wrote a script in R to plot averaged results per pairwise site distance.

```bash
# Run via UFRC queue, see plot_r2_maf0.025_long.job for more details.
# Resources: 1.54 Gb, 11 min

module load R/4.1

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/maf0.025/long_range"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for NAME in "${CHROMOSOMES[@]}"
do
        Rscript "$SCRIPTS_DIR"/plot_r2.r "$WDIR"/"$NAME"_r2.csv "$WDIR"/r2_plots/"$NAME"_r2.png "$NAME" 500
done
```

### Short Range
**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.025_short.job for more information.
# Resources: 3.95 Gb, 13 min

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.025 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 1000 --ld-window-kb 100 -ld-window-r2 0 --make-bed  --vcf-half-call m --out "$CHR" --threads 12
done
```

Commands for averaging R2 values and plotting were the same as for long range.

### Results

| Chromosome | R2 = 0.20 Distance (bp) |
| ---------- | ----------------------- |
| 1          | 265 - 399               |
| 2          | 138 - 170               |
| 3          | 72 - 89                 |
| 4          | 313 - 603               |
| 5          | 40 - 48                 |
| 6          | 307 - 474               |
| 7          | 95 - 126                |
| 8          |  105 - 122              |
| 9          |  162 - 273              |
| 10         |  318 - 776              |
| 11         |  175 - 315              |

## MAF = 0.05

### Long Range

**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.05_long.job for more information.
# Resources: 6.65 Gb, 20 min

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 100000 --ld-window-kb 1000 -ld-window-r2 0 --make-bed  --vcf-half-call m --thin 0.5 --out "$CHR" --threads 12
done
```

Commands for averaging R2 values and plotting were the same as for MAF = 0.025 long range analysis.

### Short Range
**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.05_short.job for more information.
# Resources: 3.20 Gb, 12 min

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.05 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 1000 --ld-window-kb 100 -ld-window-r2 0 --make-bed  --vcf-half-call m --out "$CHR" --threads 12
done
```

Commands for averaging R2 values and plotting were the same as for MAF = 0.025 long range analysis.

### Results
| Chromosome | R2 = 0.20 Distance (bp) |
| ---------- | ----------------------- |
| 1          | 1642 - 3969             |
| 2          | 1069 - 2603             |
| 3          | 489 - 835               |
| 4          | 2304 - 5850             |
| 5          | 317 - 567               |
| 6          | 2263 - 99723            |
| 7          | 695 - 1577              |
| 8          | 867 - 1688              |
| 9          | 1448 - 3645             |
| 10         | 2890 - 7771             |
| 11         | 1305 - 2772             |