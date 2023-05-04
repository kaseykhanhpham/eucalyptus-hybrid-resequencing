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

### "Parallelize" Chromosome 3 and 5 SNP Calling:


Chromosomes 3 and 5 required prohibitive computational resources when run in whole, so I split them up into ten parts and ran each as a separate job. Made a list of segments in which to split each chromosome.

Chromosome 3 length: 65547241
Chromosome 5 length: 62919971

Split BAM files for each sample by chromosome and segment:

```bash
# Run via job on UFRC, see split_bam.job for details
# Resources used: 12Mb, 2 hrs
module load samtools/1.15

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
CHR03_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/parallel/Chr03"
CHR05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/parallel/Chr05"

# Chromosome 3
while read NAME
do
    while read SECTION
    do
        samtools view -b -h -o "$CHR03_DIR"/bamfiles/"$NAME"_Chr03_"$SECTION".bam "$IN_DIR"/"$NAME"_marked.bam Chr03:"$SECTION"
    done < "$CHR03_DIR"/chr03_section_list.txt
done < "$LIST_DIR"/seq_ids.txt

while read SECTION
do
    samtools view -b -h -o "$CHR03_DIR"/bamfiles/SRR10339635_Chr03_"$SECTION".bam "$IN_DIR"/SRR10339635_marked.bam Chr03:"$SECTION"
done < "$CHR03_DIR"/chr03_section_list.txt

# Chromosome 5
while read NAME
do
    while read SECTION
    do
        samtools view -b -h -o "$CHR05_DIR"/bamfiles/"$NAME"_Chr05_"$SECTION".bam "$IN_DIR"/"$NAME"_marked.bam Chr05:"$SECTION"
    done < "$CHR05_DIR"/chr05_section_list.txt
done < "$LIST_DIR"/seq_ids.txt

while read SECTION
do
    samtools view -b -h -o "$CHR05_DIR"/bamfiles/SRR10339635_Chr05_"$SECTION".bam "$IN_DIR"/SRR10339635_marked.bam Chr05:"$SECTION"
done < "$CHR05_DIR"/chr05_section_list.txt
```

Manually generated separate file input lists for each chromosome section to pass to `FreeBayes`; see the files with the name template `bam_inputs_chr0X_XX.txt`.

Ran `FreeBayes` on each chromosome section separately. An example job:

```bash
# Run via job on UFRC, see fb_chr03_1.job for details
# Max resources used: 3 Gb, 6 days

module load freebayes/1.3.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/parallel/Chr03"
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/parallel/Chr03"

# Run FreeBayes on all samples at once against E. globulus reference genome, call sites with a maximum depth of 3000
# Chromosome 3, Section 1: 1 bp to 6554724 bp
freebayes -L "$LIST_DIR"/bam_inputs_chr03_01.txt -f "$REF_FILE" -v "$OUTDIR"/chr03_01.vcf -g 3000
```

Merged VCF files for each Chromosome 3/5 piece using [picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard). Example job given below.

```bash
# Run via job on UFRC, see merge_vcfs_chr03.job for details
# Max resources used: 10 Gb, 20 min

module load picard/2.25.5
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes/parallel/Chr03"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"

picard MergeVcfs -I "$VCF_DIR"/vcfs_list.txt -O "$OUT_DIR"/chr03.vcf
```

### SNP Calling Post-processing

Merged VCF files of unfiltered variants for each chromosome using `picard`.

Then calculateed raw variant statistics for the merged file:

```bash
# Run via job on UFRC, see snp_stats_raw.job for details
# Max resources used: 
module load bcftools
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"

# Get statistics for raw snps on whole genome
bcftools stats "$VCF_DIR"/meehan_all_raw.vcf > "$VCF_DIR"/raw_stats/meehan_all_raw_stats.txt
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
# Resources used: ~900 Mb, ~2 days

module load vcflib/1.0.1
module load vcftools/0.1.16
module load bcftools/1.15

export INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.freebayes"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

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
# Resources used: ~10 Mb, ~25 min

module load picard/2.25.5
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

picard MergeVcfs -I "$LIST_DIR"/vcfs_fil_to_merge.txt -O meehan_all_fil_maf0.05.vcf

# Separate SNPs and indels and remove fields that no longer fit the data:
# Run via job on UFRC, see split_indels_maf0.05.job for details
# Resources used: ~10 Gb, ~20 min

module load vcftools/0.1.16
NAME="meehan_all_fil_maf0.05"

bcftools annotate -x FORMAT/PL,FORMAT/AD,FORMAT/AO,FORMAT/QA -O v -o - "$NAME".vcf | vcftools --vcf - --remove-indels --recode --stdout > "$NAME"_snps.vcf
bcftools annotate -x FORMAT/PL,FORMAT/AD,FORMAT/AO,FORMAT/QA -O v -o - "$NAME".vcf | vcftools --vcf - --keep-only-indels --recode --stdout > "$NAME"_indels.vcf
```

Identify duplicate calls for the same position and check by eye that removed duplicates are similar to first entries for each position.

```bash
# Merge VCFs:
# Run via job on UFRC, see validate_vcf_maf0.00.job for details
# Resources used: 1 Gb, 2 min
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"

python "$SCRIPT_DIR"/check_vcf_duplicates.py "$INFILE" "meehan_all_fil_maf0.00_snps_first.txt" "meehan_all_fil_maf0.00_snps_dupes.txt" "meehan_all_fil_maf0.00_snps_nodupes.vcf"
```

## Get summary statistics for filtered variant sets

Summary statistics were generated for filtered sets in the [same way they were calculated for the raw variant set](#calculate-raw-variant-statistics). 

Variant counts after filtering:
| MAF  | SNP count | Indel count |
| ---- | --------- | ----------- |
| 0.00 | 10390327  | 2034126     |
| 0.05 | 5473125   | 984256      |

The SNP count at MAF=0.00 is much higher than without an outgroup, likely because of the inclusion of the outgroup. At MAF=0.05, the SNP counts for sets with and without outgroups are comparable.