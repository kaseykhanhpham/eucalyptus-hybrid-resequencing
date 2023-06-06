# Linkage Disequilibrium

Calculated linkage disequilibrium using r^2 decay (drop-off in correlation between SNPs of a certain distance within the genome) using [`plink v1.9`](https://www.cog-genomics.org/plink/1.9) and custom scripts.

Created a directory of symlinks to gzipped VCF files for each chromosome with the ".gz" extension removed since PLINK checks for extension propriety.

## _E. globulus_

### Calculate r^2

Calculated for all _E. globulus_ samples (both Meehan Range and reference). First, calculated r^2 values for all pairwise SNPs on each chromosome in `plink` for two different levels of minor allele frequency, since that affects LD estimations.

```bash
# Performed on UFRC queue system; see plink_glob_maf0X.job. Example from plink_glob_maf00.job given below.
# Maximum resources used: 750 Mb, 2 min

module load plink/1.90b3.39 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/vcf_files"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"
declare -a CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHROMOSOMES[@]}"
do
    plink --vcf "$IN_DIR"/"$CHR"_snps_fil.vcf --keep "$LIST_DIR"/Eglobulus.fam --double-id --allow-extra-chr --set-missing-var-ids @:# --maf 0.00000001 --mind 0.5 --r2 gz --ld-window 100001 --ld-window-kb 1000 -ld-window-r2 0 --make-bed  --vcf-half-call m --thin 0.5 --out "$CHR"_glob00 --threads 12
done
```

Averaged r^2 values for every distance between two SNPs using a custom `python` script.

```bash
# Performed on UFRC queue system; see avg_r2_glob_maf0X.job. Example from avg_r2_glob_maf00.job given below.
# Maximum resources used: 3 Gb, 30 min

module load python/3.8 

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob/maf0.00"
declare -a CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${CHROMOSOMES[@]}"
do
	python "$SCRIPTS_DIR"/average_r2.py  "$WDIR"/"$NAME"_glob00.ld.gz "$WDIR"/"$NAME"_r2_glob00.csv
done
```

Averaged r^2 values across chromosomes using a custom `python` script.

```bash
# Performed on UFRC queue system; see genwide_glob_avg_maf0X.job for more details. Example from genwide_glob_avg_maf00.job given below.
# Maximum resources used: 250 Mb, 7 hrs

module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts" 

python "$SCRIPT_DIR"/genomewide_r2_avg.py {CHR}_r2_glob00.csv
```

### Process r^2 calculations

Plotted r^2 using a custom `R` script.

```bash
# Performed on UFRC queue system; see plot_r2_glob_maf0X.job for more details. Example from plot_r2_glob_maf00.job given below.
# Maximum resources used: 1.2 Gb, 30 min

module load R/4.1 

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
declare -a CHROMOSOMES=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for NAME in "${CHROMOSOMES[@]}"
do
	Rscript "$SCRIPTS_DIR"/plot_r2.r "$NAME"_r2_glob00.csv r2_plots/"$NAME"_r2_glob00.png "$NAME" 500
done
```

Fitted r^2 decay curve. ???? IDK how...

## _E. cordata_

Followed the same steps for _E. cordata_ samples. See jobfiles with the "cord" tag for more details.