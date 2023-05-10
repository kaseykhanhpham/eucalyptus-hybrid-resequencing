# Genome Scan Analysis

## Calculate statistics
Calculated diagnostic statistics for introgression in sliding windows across genome using the python library `egglib`. Wrote a script to automate egglib calculations for sliding windows: `stats_windows.py`. It requires the following input files present in this directory:

`chr_list.txt`: list of chromosome/contig names in reference genome, one per line

**structure files**: The JSON files in this directory define the individuals and populations to be used for calculations of diagnostic stats. The "glob" files consider introgressed and reference _E. globulus_ to be the ingroup (stats calculated) and _E. cordata_ and _E. grandis_ to be the outgroup (stats not calculated). The "cord" files consider introgressed _E. globulus_ and _E. cordata_ to be the ingroup and reference _E. globulus_ and _E. grandis_ to be the outgroup.

Ran `stats_windows.py` with the stats arguments "Pi", "FST", "Dxy", "Deta" (see below).

```bash
# Done on UFRC queue system, see pi_windows.job, fst_windows.job, dxy_windows.job, deta_windows.job for more detail
# Code below is example from pi_windows.job
# Resources used:

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/stats_windows.py "$VCF_DIR"//meehan_all_fil_maf0.00_snps.vcf Pi 100000 20000 "$WDIR"/glob_structure.json "$WDIR"/glob_output.json "$WDIR"/chr_list.txt glob_pi_windows.tab
```

## Patterson's D (and associated)

First filtered VCF to only biallelic or monoallelic sites.

```bash
# Performed in UFRC queue system. See get_biallelic.job for more details.
# Resources: 3 Mb, 13 min
module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00"

vcftools --vcf "$VCFDIR"/meehan_all_fil_maf0.00_snps.vcf --min-alleles 2 --max-alleles 2 --recode --stdout > meehan_all_fil_maf0.00_snps_biallelic.vcf
```

Used biallelic SNP set to calculate D statistics for individuals and in a sliding window across genomes, using _E. grandis_ as the outgroup.

```bash
# Performed in UFRC queue system. See dsuite.job for more details.
# Resources: 4 Mb, 17 min
DSUITE_DIR="/blue/soltis/kasey.pham/bin/Dsuite/Build"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"
NAME="meehan_all_fil_maf0.00_snps_biallelic"

"$DSUITE_DIR"/Dsuite Dinvestigate -n "$NAME" "$WDIR"/"$NAME".vcf "$WDIR"/SETS.txt test_trios.txt
```
