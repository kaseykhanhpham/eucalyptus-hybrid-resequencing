# Linkage Disequilibrium Calculation

DON'T DO THIS IMMEDIATELY... REVISE TO FIT DECAY CURVE FOR MORE FORMALIZED ESTIMATE.
[Meier and Ravinet's 2019 Speciation Genomics tutorial](https://speciationgenomics.github.io/ld_decay/) referenced extensively. LD only calculated on main chromosomes because smaller contigs did not contain enough SNPs to be informative.

## LD for E. globulus

### Calculate pairwise variance between SNPs

```bash
# Run on UFRC's queue system, see plink_maf0.0001.job and plink_maf0.05.job for more information.
# plink_maf0.0001.job example commands below; only maf was changed between jobs.
# Max resource usage: 22 Gb, 50 min

module load plink/1.90b3.39 

export INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$INFILE" --keep "$LIST_DIR"/Eglobulus.fam --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.00001 --mind 0.5 --chr "$CHR" --r2 gz --ld-window 100001 --ld-window-kb 1000 -ld-window-r2 0 --make-bed  --vcf-half-call m --thin 0.5 --out "$CHR" --threads 12
done
```

**Process plink outputs:**

Wrote python script to average pairwise R2 values for each distance between SNPs and ran for each chromosome.
```bash
# Run via UFRC queue, see avg_r2_maf0.0001.job and avg_r2_maf0.05 for more details.
# avg_r2_maf0.00001.job shown below
# Max resources used: 110 Gb, 15 hrs

module load python/3.8 

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/maf0.00001"
CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${CHROMOSOMES[@]}"
do
	python "$SCRIPTS_DIR"/average_r2.py  "$WDIR"/"$NAME".ld.gz "$WDIR"/"$NAME"_r2.csv
done
```

Average across all chromosomes in `python`:

```bash
# Performed in UFRC queue system; see genomewide_avg_maf0.00001.job and genomewide_avg_maf0.05 for more details
# Example before
# Resources used: 220 Mb, 4 hrs
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts" 

python "$SCRIPT_DIR"/genomewide_r2_avg.py
```

**Plot plink output:**

Wrote a script in R to plot averaged results per pairwise site distance.

```bash
module load R/4.1

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/maf0.025/long_range"
CHROMOSOMES=(NC_052612.1 NC_052613.1 NC_052614.1 NC_052615.1 NC_052616.1 NC_052617.1 NC_052618.1 NC_052619.1 NC_052620.1 NC_052621.1 NC_052622.1)

for NAME in "${CHROMOSOMES[@]}"
do
        Rscript "$SCRIPTS_DIR"/plot_r2.r "$WDIR"/"$NAME"_r2.csv "$WDIR"/r2_plots/"$NAME"_r2.png "$NAME" 500
done
```

Genome-wide average results:

!["Linkage Disequilibrium MAF0.00001, average r2 between positions vs. distance between loci (bp)"](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/ld/glob/maf0.00001/genomewide_r2_plot_full.png "Linkage Disequilibrium MAF0.00001, Full")

!["Linkage Disequilibrium MAF0.00001, average r2 between positions vs. distance between loci (bp)"](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/ld/glob/maf0.00001/genomewide_r2_plot_zoomed.png "Linkage Disequilibrium MAF0.00001, Zoomed")
