# Linkage Disequilibrium Calculation

DON'T DO THIS IMMEDIATELY... REVISE TO FIT DECAY CURVE FOR MORE FORMALIZED ESTIMATE.
[Meier and Ravinet's 2019 Speciation Genomics tutorial](https://speciationgenomics.github.io/ld_decay/) referenced extensively. LD only calculated on main chromosomes because smaller contigs did not contain enough SNPs to be informative.

## MAF = 0.025
### Long Range

**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.025_long.job for more information.
# Resources: 4 to 13 Gb, 10 to 40 min

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.025 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 100000 --ld-window-kb 1000 -ld-window-r2 0 --make-bed  --vcf-half-call m --thin 0.5 --out "$CHR" --threads 12
done
```

**Process plink outputs:**

Wrote python script to average pairwise R2 values for each distance between SNPs and ran for each chromosome.
```bash
# Run via UFRC queue, see avg_r2_maf0.025_long.job for more details.
# Resources: 21 Gb to 65 Gb, 3 to 10 hrs

module load python/3.8

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/maf0.025/long_range"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${CHROMOSOMES[@]}"
do
        python "$SCRIPTS_DIR"/average_r2.py  "$WDIR"/"$NAME".ld.gz "$WDIR"/"$NAME"_r2.csv
done
```

**Plot plink output:**

Wrote a script in R to plot averaged results per pairwise site distance.

```bash
# Run via UFRC queue, see plot_r2_maf0.025_long.job for more details.
# Resources: 300 Mb to 1.5 Gb, 2 to 20 min

module load R/4.1

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/maf0.025/long_range"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for NAME in "${CHROMOSOMES[@]}"
do
        Rscript "$SCRIPTS_DIR"/plot_r2.r "$WDIR"/"$NAME"_r2.csv "$WDIR"/r2_plots/"$NAME"_r2.png "$NAME" 500
done
```

The same was done to calculate short-range LD using a window of 1000bp and a maximum distance of 100kb. 

### Short Range
**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink_maf0.025_short.job for more information.

module load plink/1.90b3.39

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.025 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 1000 --ld-window-kb 100 -ld-window-r2 0 --make-bed  --vcf-half-call m --out "$CHR" --threads 12
done
```

Commands for averaging R2 values and plotting were the same as for long range.

**Summarize LD results using custom python script:**
```bash
module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHROMOSOMES[@]}"
do
    echo doing "$CHR" >> avg_r2_summary.txt
    python "$SCRIPTS_DIR"/summarize_r2.py "$CHR"_r2.csv 0.2 >> avg_r2_summary.txt
done
```

### Results for MAF = 0.025 at a cutoff of r2 = 0.2

**Long Range:**
| Chromosome | Minimum (bp) | Maximum (bp) | Midpoint (bp) |
| ---------- | ------------ | ------------ | ------------- |
| 1          | 268          | 670          | 469           |
| 2          | 162          | 236          | 199           |
| 3          | 81           | 107          | 94            |
| 4          | 384          | 831          | 607.5         |
| 5          | 44           | 61           | 52.5          |
| 6          | 362          | 733          | 547.5         |
| 7          | 105          | 149          | 127           |
| 8          |  120         | 136          | 128           |
| 9          |  194         | 360          | 277           |
| 10         |  405         | 1369         | 887           |
| 11         |  228         | 433          | 330.5         |
| Average    |  -           |    -         | 338.1         |

## MAF = 0.05
Repeated the above pipeline for a cutoff of MAF=0.05 for long and short range windows. Commands for calculating pairwise LD, averaging over SNP distance, and summarizing were the same as described above. See individual job files for details.

### Results
| Chromosome | Minimum (bp) | Maximum (bp) | Midpoint (bp) |
| ---------- | ------------ | ------------ | ------------- |
| 1          | 2310         | 5614         | 3962          |
| 2          | 1272         | 3740         | 2506          |
| 3          | 599          | 1657         | 1128          |
| 4          | 2733         | 7357         | 5045          |
| 5          | 375          | 647          | 511           |
| 6          | 2831         | 6662         | 4746.5        |
| 7          | 799          | 2119         | 1459          |
| 8          | 954          | 2137         | 1545.5        |
| 9          | 2056         | 5438         | 3747          |
| 10         | 3594         | 9001         | 6297.5        |
| 11         | 1622         | 3663         | 2642.5        |
| Average    |  -           | -            | 3053.6        |