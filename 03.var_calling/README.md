# Variant Calling

## Call variants from mapped reads

Called variants in [`bcftools mpileup`](https://samtools.github.io/bcftools/howtos/variant-calling.html) to generate an AllSites VCF file. Called each chromosome in tenths to merge back together afterwards.

```bash
# Done in UFRC queue system. See mp_chrXX_X.job for more details. Example from mp_chr01_1.job below.
# Maximum resources used: 250 Mb, 1 hr

module load bcftools/1.15

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup/chr01"

bcftools mpileup --bam-list "$LIST_DIR"/bam_inputs.txt --full-BAQ --max-depth 3000 --fasta-ref "$REF_FILE" --regions Chr01:1-4221956 --seed 9387492 --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP --output-type u --threads 12 | bcftools call --ploidy 2 --output-type z --output "$OUTDIR"/chr01_1.vcf.gz --threads 12 --annotate FORMAT/GQ --multiallelic-caller
```

Merge VCF chunks for each chromosome using `bcftools`:

```bash
# Done in UFRC queue system. See mp_chrXX_X.job for more details. Example from merge_vcfs_chr01.job below.
# Maximum resources used: 25 Mb, 20 min

module load bcftools/1.15
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup/chr01"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"

ls "$VCF_DIR"/*.vcf.gz > "$VCF_DIR"/vcfs_list.txt
bcftools concat -f "$VCF_DIR"/vcfs_list.txt -O z -o "$OUT_DIR"/chr01.vcf.gz --threads 12
```

Get stats for raw variant sets using `bcftools`:

```bash
# Done in UFRC queue system. See raw_stats.job for more details.
# Resources used: 8 Mb, 35 min

module load bcftools/1.15
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup/raw_stats"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    bcftools stats --fasta-ref "$REF_FILE" "$VCF_DIR"/"$NAME".vcf.gz > "$OUT_DIR"/"$NAME"_raw_stats.txt
done
```

Visualize raw variant set stats summary using script from `bcftools`:
```bash
# Done in UFRC queue system. See vis_stats.job for more details.
# Resources used: 51 Mb, 1 min
module load bcftools/1.15
module load python/3.10
module load texlive

REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup/raw_stats"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    "$HPC_BCFTOOLS_DIR"/bin/plot-vcfstats -p "$OUT_DIR"/"$NAME"_vis -t "$NAME" -v "$NAME"_raw_stats.txt
done
```

## Filter variants

Split raw variant calls into invariant SNPs, variant SNPs, and indels using `vcftools` to assess which filtering parameters to use. Ignoring indels in analysis for now to focus on SNPs.

```bash
# Done in UFRC queue system. See split_vcfs.job
# Maximum resources used:
module load vcftools/0.1.16
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --remove-indels --max-maf 0 --recode --stdout | bgzip -c > "$NAME"_snp_invar.vcf.gz
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --remove-indels --mac 1 --recode --stdout | bgzip -c > "$NAME"_snp_var.vcf.gz
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --keep-only-indels --recode --stdout | bgzip -c > "$NAME"_indels.vcf.gz
done
```

