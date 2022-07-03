# Linkage Disequilibrium Calculation

[Meier and Ravinet's tutorial](https://speciationgenomics.github.io/ld_decay/) referenced extensively.

**Run plink to get pairwise r2 values:**
```bash
# Run on UFRC's queue system, see plink.job for more information.
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

