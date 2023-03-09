# Admixture Mapping
Admixture mapping is the process of assigning loci within the genome of an admixed individual to parental populations using allele frequencies in a phased dataset.

Overview of the method: [Winkler et al. 2010](https://doi.org/10.1146/annurev-genom-082509-141523)


## Phase Alleles
### Remove variants with missing genotypes
```bash
# Performed in UFRC queue system. See filter_missing.job for more details.
# Resources: 75 Mb, 1 min

module load vcftools/0.1.16
module load bcftools/1.15

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05"
NAME="all_to_ASM1654582_fil_maf0.05_snps"

export MISS=1.0

cat "$INDIR"/"$NAME".vcf | vcftools --vcf - --max-missing $MISS --recode --stdout | bcftools sort -O v - > "$NAME"_nomiss.vcf
```

### Check VCF for invalid entries using [Xiaowei Zhan's checkVCF.py script](https://github.com/zhanxw/checkVCF/blob/master/checkVCF.py)
No indications from check programs what the issue was with my VCF file here.

```bash
module load gatk
gatk ValidateVariants -R /blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1/GCF_016545825.1_ASM1654582v1_genomic.fna -V all_to_ASM1654582_fil_maf0.05_snps_nomiss.vcf > validate_vcf.txt # no issues

module load python/2

SCRIPTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
GENOMEDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"

python "$SCRIPTDIR"/checkVCF.py -r "$GENOMEDIR"/GCF_016545825.1_ASM1654582v1_genomic.fna -o all_to_ASM1654582_fil_maf0.05_snps_nomiss all_to_ASM1654582_fil_maf0.05_snps_nomiss.vcf
```

### Filter all variants with samples with haploid GT calls
```bash
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
module load bcftools/1.15

bcftools filter -e 'FORMAT/GT="hap"' -O v -o all_to_ASM1654582_fil_maf0.05_snps_nohap.vcf all_to_ASM1654582_fil_maf0.05_snps_nomiss.vcf
```

### Phase alleles with [`Beagle`](https://faculty.washington.edu/browning/beagle/beagle.html)

```bash
# Performed in UFRC queue system. See filter_missing.job for more details.
# Resources: 1.52 Mb, 30s
module load beagle/5.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/05.analyses/admix_map"
NAME="all_to_ASM1654582_fil_maf0.05_snps_nohap"

"$HPC_BEAGLE_BIN"/beagle gt="$NAME".vcf out="$NAME"_phased impute=false burnin=5 iterations=20 phase-states=280 nthreads=11
```

### Run admixture mapping software [`flare`](https://github.com/browning-lab/flare)

```bash

```