# Admixture Mapping
Admixture mapping is the process of assigning loci within the genome of an admixed individual to parental populations using allele frequencies in a phased dataset.

Overview of the method: [Winkler et al. 2010](https://doi.org/10.1146/annurev-genom-082509-141523)


## Phase Alleles
### Format VCF for phasing in BEAGLE
Removed variants with missing genotypes, invariants, genotypes with haploid calls, and unanchored contigs. Restricted to only ingroup samples.

```bash
# Run on UFRC queue system; see format_beagle.job for more details.
# Resources used: 2 Mb, 1 min

module load vcftools/0.1.16
module load bcftools/1.15

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
NAME="all_fil"

export MISS=1.0

vcftools --gzvcf "$INDIR"/"$NAME".vcf.gz --max-missing $MISS --min-alleles 2 --not-chr ChrUn --recode --stdout | bcftools view -S "$LIST_DIR"/sample_ids.txt - | bcftools filter -e 'FORMAT/GT="hap"' -O v - | bcftools sort -O v - > "$NAME"_beagle_formatted.vcf
```

Phased variants.

```bash
# Run on UFRC queue system; see beagle.job for more details.
# Resources used: 2 Mb, 10 sec

module load beagle/5.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/05.analyses/admix_map"
NAME="all_fil"

"$HPC_BEAGLE_BIN"/beagle gt="$NAME"_beagle_formatted.vcf out="$NAME"_phased impute=false burnin=5 iterations=20 phase-states=280 nthreads=11
```

Ran `flare` to do admixture mapping.

```bash
# Run on UFRC queue system; see flare.job for more details.
# Resources used: 2 Mb, 5 sec

module load java/11.0.1

FLARE_DIR="/blue/soltis/kasey.pham/bin/flare"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
NAME="all_fil"

java -Xmx10g -jar "$FLARE_DIR"/flare.jar ref="$WDIR"/"$NAME"_phased.vcf.gz ref-panel="$WDIR"/ref_panel.txt gt="$WDIR"/"$NAME"_phased.vcf.gz gt-samples="$WDIR"/gt_samples.txt map="$WDIR"/1060_LH_F2_manual_copy.map min-mac=1 nthreads=12 out="$NAME"_flare
```
