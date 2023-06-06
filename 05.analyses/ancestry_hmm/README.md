# Local Ancestry Inference
# HMM Inference of Local Ancestry
Performed using the tool [`Ancestry_HMM`](https://github.com/russcd/Ancestry_HMM). Thanks to Shelley Sianta for suggesting this analysis.

## Prepare inputs
### Genotype count input file

First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts.job for more details.
# Resources used: 5 Mb, 50 min

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out "$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in.job for more details.
# Resources used: 5 Mb, 50 min
```