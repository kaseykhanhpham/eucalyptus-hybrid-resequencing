# Admixture Mapping
Admixture mapping is the process of assigning loci within the genome of an admixed individual to parental populations using allele frequencies in a phased dataset.

Overview of the method: [Winkler et al. 2010](https://doi.org/10.1146/annurev-genom-082509-141523)



## Phase Alleles
### Remove variants with missing genotypes
```bash
# Performed in UFRC queue system. See filter_missing.job for more details.
# Resources: 80 Mb, 2 min

module load vcftools/0.1.16
module load bcftools/1.15

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05"
NAME="all_to_ASM1654582_fil_maf0.05"

export MISS=1.0

cat "$INDIR"/"$NAME".vcf | vcftools --vcf - --max-missing $MISS --recode --stdout | bcftools sort -O v - > "$NAME"_nomiss.vcf
```

### Check VCF for invalid entries using [Xiaowei Zhan's checkVCF.py script](https://github.com/zhanxw/checkVCF/blob/master/checkVCF.py)


```bash
module load python/2

SCRIPTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
GENOMEDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"

python "$SCRIPTDIR"/checkVCF.py -r "$GENOMEDIR"/GCF_016545825.1_ASM1654582v1_genomic.fna -o all_to_ASM1654582_fil_maf0.05_nomiss all_to_ASM1654582_fil_maf0.05_nomiss.vcf
```

### Phase alleles with [`Beagle`](https://faculty.washington.edu/browning/beagle/beagle.html)

```bash
# Performed in UFRC queue system. See filter_missing.job for more details.
# Resources: 

```