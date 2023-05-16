# Variant Calling

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

Merge VCF chunks for each chromosome:

```bash
# Done in UFRC queue system. See mp_chrXX_X.job for more details. Example from merge_vcfs_chr01.job below.
# Maximum resources used: 25 Mb, 20 min

module load bcftools/1.15
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup/chr01"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"

ls "$VCF_DIR"/*.vcf.gz > "$VCF_DIR"/vcfs_list.txt
bcftools concat -f "$VCF_DIR"/vcfs_list.txt -O z -o "$OUT_DIR"/chr01.vcf.gz --threads 12
```
