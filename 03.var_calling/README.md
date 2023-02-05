# Variant Calling

## Housekeeping:
```bash
READ_INDIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/04.processed_reads"
INDEX_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
READ_OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
VC_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"

# create symlinks to processed reads in variant calling directory
while read NAME; do ln -s "$READ_INDIR"/"$NAME"_marked.bam "$READ_OUTDIR"/"$NAME"_marked.bam; ln -s "$INDEX_DIR"/"$NAME"_marked.bai "$READ_OUTDIR"/"$NAME"_marked.bai; done < "$LIST_DIR"/seq_ids.txt
# create list of BAM inputs for variant caller
ls -d "$READ_OUTDIR"/* > "$VC_DIR"/bam_inputs.txt
```

## Run Variant Caller
I am using the variant caller [`FreeBayes`](https://github.com/freebayes/freebayes), which estimates internal population parameters for the individuals it is provided with using Bayesian inference. This is probably a good choice for my dataset given that I don't have a set of validated SNPs for training a variant caller like `GATK`. I followed recommendations on `FreeBayes`' Github page for parameter selection.

I ran `FreeBayes` on each chromosome separately, with the same parameters each time. See the job files named `fb_chr--.job` for each set of individual parameters.

**Command for chromosome 1:**

```bash
# Run via job on UFRC, see fb_chr1.job for details
# Resources used: 14 days, 5.35 Gb
# (range: 4 Gb - 8 Gb, 11 days - 27 days)
# Some chromosomes did not run until I gave them 20 Gb of RAM

module load freebayes/1.3.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/"

# Run FreeBayes on all samples at once against E. globulus reference genome, call sites with a maximum depth of 4000
freebayes -L "$LIST_DIR"/bam_inputs.txt -f "$REF_FILE" -r Chr01 -v "$OUTDIR"/Chr01.vcf -g 4000
```

**Merge VCF files for each chromosome using [picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard):**

```bash
# Run via job on UFRC, see merge_vcfs_raw.job for details
# Resources used: 10.2 Gb, 2 hrs 

module load picard/2.25.5
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"

ls "$VCF_DIR"/*.vcf > "$VCF_DIR"/vcfs_list.txt
picard MergeVcfs -I "$VCF_DIR"/vcfs_list.txt -O meehan_all_raw.vcf
```

Then calculateed raw variant statistics for the merged file:

```bash
module load bcftools
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"

# Get statistics for entire variant set
bcftools stats "$VCF_DIR"/meehan_all_raw.vcf > "$VCF_DIR"/all_raw_stats.txt
```

## Filter SNPs
SNP filtering parameters and code were based on several sources, including personal consultation with Zhe Cai (University of British Columbia), Robin Buell's Fall 2018 Plant Genomics course at Michigan State University, [Hübner et al. 2019](https://doi.org/10.1038/s41477-018-0329-0), and Ravinet and Meier's [Speciation and Population Genomics](https://speciationgenomics.github.io/) tutorial.

**Filtering Parameters used:**

| Filter / Command                            | Meaning                                                             |
| ------------------------------------------- | ------------------------------------------------------------------- |
| [vcfallelicprimitives](https://github.com/vcflib/vcflib/blob/master/doc/vcfallelicprimitives.md) | Separate multi-nucleotide polymorphisms (MNPs) into SNPs |
| SAF > 1 & SAR > 1                           | Remove variants with strand bias                                    |
| RPR > 1 & RPL > 1                           | Remove variants with position bias relative to reference allele     |
| maf = 0.05                                  | Remove variants with a minor allele frequency < 0.05                |
| max-missing = 0.9                           | Remove variants with a more than 10% missing genotypes              |
| minQ = 50                                   | Remove variants with a quality (QUAL) score < 50                    |
| min-meanDP = 10,                            | Remove variants with an average per-sample depth (DP) < 10          |
| max-meanDP = 40                             | Remove variants with an average per-sample depth (DP) > 40          |
| minDP = 10                                  | Remove sample genotypes with a depth (DP) < 10                      |
| maxDP = 40                                  | Remove sample genotypes with a depth (DP) > 40                      |
| recode                                      | Write header in VCF format                                          |
| QUAL/MNT/DP < 20                            | Remove variants with quality normalized by depth (QUAL/DP) < 20     |
| FMT/AQ/FMT/AO < 30                          | Remove variants with average alternate allele quality (AQ/AO) < 30  |
| sort                                        | Sort filtered variants by position                                  |

**Programs used:**

* [`vcflib`](https://github.com/vcflib/vcflib)
* [`vcftools`](https://vcftools.github.io)
* [`bcftools`](https://samtools.github.io/bcftools/bcftools.html)


Four filtering sets were generated in total, only varying in minimum threshold for minor allele frequencing (MAF): {MAF > 0.00, 0.05, 0.10, 0.30}. Command for `MAF > 0.05` dataset shown below, but see job files `filtersnps_maf---.job` and `merge_vcfs_maf---.job` for more detail.

**Command for `MAF > 0.05` and merging of filtered chromosome VCFs:**

```bash
# Filter SNPs:
# Run via job on UFRC, see filtersnps_maf0.05.job for details
# Resources used: 900 Mb, 2 days (range: 400 Mb - 1 Gb, 35 hrs - 48 hrs)

module load vcflib/1.0.1
module load vcftools/0.1.16
module load bcftools/1.15

export INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes/max_dp_3000"
declare -a VCFLIST=(NC_052612 NC_052613 NC_052614 NC_052615 NC_052616 NC_052617 NC_052618 NC_052619 NC_052620 NC_052621 NC_052622 unanchored_contigs)

export MAF=0.05
export MISS=0.9
export QUAL=50
export MIN_DP=10
export MAX_DP=40

for NAME in "${VCFLIST[@]}"
do
    echo doing "$NAME"
    cat "$INDIR"/"$NAME".vcf | vcfallelicprimitives -kg | vcffilter -f "SAF > 1 & SAR > 1 & RPR > 1 & RPL > 1" | vcftools --vcf - --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DP --max-meanDP $MAX_DP --minDP $MIN_DP --maxDP $MAX_DP --recode --stdout | bcftools view -e 'QUAL/FMT/DP<20 & FMT/QA/FMT/AO<30' -O v - | bcftools sort -O v - > "$NAME"_fil.vcf
done

# Merge VCFs:
# Run via job on UFRC, see merge_vcfs_maf0.05.job for details
# Resources used: 10.2 Gb, 10 min (range: 10.2 Gb, 2 - 10 min)

module load picard/2.25.5
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps"

# Create list of VCF files to merge
ls *.vcf > "$LIST_DIR"/vcfs_fil.list
# Merge VCFs
picard MergeVcfs -I "$LIST_DIR"/vcfs_fil.list -O all_to_ASM1654582_fil.vcf

# Separate SNPs and indels:
# Run via job on UFRC, see split_indels_maf0.05.job for details
# Resources used: 8.5 Mb, 25 min

module load vcftools/0.1.16
NAME="all_to_ASM1654582_fil_maf0.05"

cat "$NAME".vcf | vcftools --vcf - --remove-indels --recode --stdout > "$NAME"_snps.vcf
cat "$NAME".vcf | vcftools --vcf - --keep-only-indels --recode --stdout > "$NAME"_indels.vcf
```
## Get summary statistics for filtered variant sets

Summary statistics and visualizations were generated for filtered sets in the [same way they were calculated for the raw variant set](#calculate-and-visualize-raw-variant-statistics). See `vis_chr_maf---.job` files for more detail. 

**Calculate genome-wide statistics:**
```bash
module load bcftools
MAF00DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00"
MAF05DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05"
MAF10DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.10"
MAF30DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.30"

bcftools stats "$MAF00DIR"/all_to_ASM1654582_fil_maf0.00_snps.vcf > "$MAF00DIR"/all_fil_maf0.00_snp_stats.txt
bcftools stats "$MAF05DIR"/all_to_ASM1654582_fil_maf0.05_snps.vcf > "$MAF05DIR"/all_fil_maf0.05_snp_stats.txt
bcftools stats "$MAF10DIR"/all_to_ASM1654582_fil_maf0.10_snps.vcf > "$MAF10DIR"/all_fil_maf0.10_snp_stats.txt
bcftools stats "$MAF30DIR"/all_to_ASM1654582_fil_maf0.30_snps.vcf > "$MAF30DIR"/all_fil_maf0.30_snp_stats.txt

bcftools stats "$MAF00DIR"/all_to_ASM1654582_fil_maf0.00_indels.vcf > "$MAF00DIR"/all_fil_maf0.00_indel_stats.txt
bcftools stats "$MAF05DIR"/all_to_ASM1654582_fil_maf0.05_indels.vcf > "$MAF05DIR"/all_fil_maf0.05_indel_stats.txt
bcftools stats "$MAF10DIR"/all_to_ASM1654582_fil_maf0.10_indels.vcf > "$MAF10DIR"/all_fil_maf0.10_indel_stats.txt
bcftools stats "$MAF30DIR"/all_to_ASM1654582_fil_maf0.30_indels.vcf > "$MAF30DIR"/all_fil_maf0.30_indel_stats.txt
```

Chromosome and genome-wide statistics for raw and filtered variant sets can be found in the subdirectory [variant_statistics](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/03.var_calling/variant_statistics) and have been compiled in the file [mdp3000.xlsx](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/03.var_calling/variant_statistics/mdp3000.xlsx).
