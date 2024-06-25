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

Get stats for raw genotype call sets using `bcftools`:

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

Visualize raw genotype call set stats summary using script from `bcftools`:
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

## Filter calls
### Evaluate distributions of annotations for raw SNPs before setting filtering cutoffs

Split raw genotype calls into invariant SNPs, variant SNPs, and indels using `vcftools`. Assessed which filtering parameters to use. Ignoring indels in analysis for now to focus on SNPs.

```bash
# Done in UFRC queue system. See split_vcfs.job for more detail.
# Resources used: 15 Mb, 20 hrs

module load vcftools/0.1.16 
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/03.mpileup"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --max-maf 0 --recode --stdout | bgzip -c > "$NAME"_invar.vcf.gz
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --mac 1 --remove-indels --recode --stdout | bgzip -c > "$NAME"_snps.vcf.gz
    vcftools --gzvcf "$VCF_DIR"/"$NAME".vcf.gz --mac 1 --keep-only-indels --recode --stdout | bgzip -c > "$NAME"_indels.vcf.gz
done
```

Retrieve annotations for SNP variants.

```bash
# Done in UFRC queue system. See explore_snp_stats.job for more detail.
# Resources used: 5 Mb, 50 min
module load bcftools/1.15
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/raw_split"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

SAMPLE_HEADER="WB02\tWC02\tWD02\tWE02\tWF02\tWG02\tWH02\tWA03\tWB03\tWC03\tWA01\tWD03\tWE03\tWF03\tWG03\tWH03\tWA04\tWB04\tWC04\tWD04\tWE04\tWB01\tWF04\tWG04\tWC01\tWD01\tWE01\tWF01\tWG01\tWH01\tWA02\tWH04\tWA05\tWB05\tWC05\tWD05\tWE05\tWF05\tWG05\tWH05\n"

for NAME in "${VCFLIST[@]}"
do
    # INFO fields
    printf 'CHROM\tPOS\tQUAL\tMQ\tMQBZ\tBQBZ\n' > "$NAME"_variant_info.tab
    bcftools query --format '%CHROM\t%POS\t%QUAL\t%MQ\t%MQBZ\t%BQBZ\n' "$VCF_DIR"/"$NAME"_snps.vcf.gz >> "$NAME"_variant_info.tab
    # DP
    printf "$SAMPLE_HEADER" > "$NAME"_variant_dp.tab
    bcftools query --format '[ %DP]\n' "$VCF_DIR"/"$NAME"_snps.vcf.gz >> "$NAME"_variant_dp.tab
    # SP
    printf "$SAMPLE_HEADER" > "$NAME"_variant_sp.tab
    bcftools query --format '[ %SP]\n' "$VCF_DIR"/"$NAME"_snps.vcf.gz >> "$NAME"_variant_sp.tab
    # GQ
    printf "$SAMPLE_HEADER" > "$NAME"_variant_gq.tab
    bcftools query --format '[ %GQ]\n' "$VCF_DIR"/"$NAME"_snps.vcf.gz >> "$NAME"_variant_gq.tab
done
```

Plot various annotations in `R`:

```bash
# Done in UFRC queue system. See plot_snp_stats.job for more detail.
# Resources used: 8 Gb, 50 min
module load R/4.2

declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

for NAME in "${VCFLIST[@]}"
do
    Rscript "$SCRIPT_DIR"/plot_snp_stats.r -d "$NAME"_variant_dp.tab -g "$NAME"_variant_gq.tab -s "$NAME"_variant_sp.tab -i "$NAME"_variant_info.tab -o "$NAME"
done
```

Plotted zoomed views of stats with heavy skews.

```R
# Chr01 SP (as a proxy for the whole genome)
chr01_sp <- read.table("chr01_variant_sp.tab", header = TRUE, na.strings=".")
sp_avg <- apply(chr01_sp, 1, mean)
sp_avg_trunc <- sort(sp_avg)[which(sort(sp_avg) < 1)]
png("chr01_sp_distr_zoomed.png", width = 1000, height = 1000, unit = "px")
    hist(sp_avg_trunc, main = "Chr01 Average SP Distribution Zoomed", xlab = "Average SP")
dev.off()

# Chr03 DP (as a proxy for the whole genome)
chr03_dp <- read.table("chr03_variant_dp.tab", header = TRUE, na.strings=".")
dp_avg <- apply(chr03_dp, 1, mean)
dp_avg_trunc <- sort(dp_avg)[which(sort(dp_avg)<101)]
png("chr03_dp_distr_avg.png", width = 1000, height = 1000, unit = "px")
    hist(dp_avg, main = "Chr03 Average Depth Distribution", xlab = "Average DP")
dev.off()
png("chr03_dp_distr_avg_zoomed.png", width = 1000, height = 1000, unit = "px")
    hist(dp_avg_trunc, main = "Chr03 Average Depth Distribution Zoomed", xlab = "Average DP")
dev.off()

# Chr03 QUAL (as a proxy for the whole genome)
dp_total <- apply(chr03_dp, 1, sum)
chr03_qual <- read.table("chr03_variant_info.tab", header = TRUE, na.strings=".")$QUAL
qual_corr <- chr03_qual / dp_total
png("chr03_qual_corrected_distr.png", width = 1000, height = 1000, unit = "px")
    hist(qual_corr, main = "Chr03 Quality (Depth-Corrected) Distribution", xlab = "QUAL / DP")
dev.off()
qual_corr_trunc <- sort(qual_corr)[which(sort(qual_corr) < 2.5)]
png("chr03_qual_corrected_distr_zoomed.png", width = 1000, height = 1000, unit = "px")
    hist(qual_corr_trunc, main = "Chr03 Quality (Depth-Corrected) Distribution Zoomed", xlab = "QUAL / DP")
dev.off()
```

### Filter called SNPs

SNP filtering parameters and code were based on several sources, including personal consultation with Zhe Cai (University of British Columbia), Robin Buell's Fall 2018 Plant Genomics course at Michigan State University, [Hübner et al. 2019](https://doi.org/10.1038/s41477-018-0329-0), and Ravinet and Meier's [Speciation and Population Genomics](https://speciationgenomics.github.io/) tutorial.

**Filtering Parameters used for Variants:**

| Filter / Command                            | Meaning                                                               |
| ------------------------------------------- | --------------------------------------------------------------------- |
| maf = 0.00                                  | No filtering based on minor allele frequency at this time.            |
| minQ = 50                                   | Remove variants with a quality (QUAL) score < 50                      |
| min-meanDP = 20,                            | Remove variants with an average per-sample depth (DP) < 20            |
| max-meanDP = 60                             | Remove variants with an average per-sample depth (DP) > 60            |
| minDP = 20                                  | Remove sample genotypes with a depth (DP) < 10                        |
| maxDP = 60                                  | Remove sample genotypes with a depth (DP) > 50                        |
| recode                                      | Write header in VCF format                                            |
| QUAL/SMPL_SUM(AD) > 20                      | Remove gts with a ratio of QUAL to total allelic depth per sample < 20|
| SP > 0.10                                   | Remove genotypes with phred-corrected P-value for strand bias < 0.10  |
| GQ > 40                                     | Remove genotypes with a genotype quality < 40                         |
| max-missing = 0.875                         | Remove variants with a more than 12.5% missing genotypes              |
| sort                                        | Sort filtered variants by position                                    |

**Programs used:**

* [`vcftools`](https://vcftools.github.io)
* [`bcftools`](https://samtools.github.io/bcftools/bcftools.html)

Filtered variants first to adjust parameters to get enough SNPs to use later.

```bash
# Done via UFRC queue system; see filter_vars.job for more details.
# Resources used: 525 Mb, 40 min

module load vcftools/0.1.16
module load bcftools/1.15

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/raw_split"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

export MAF=0.00
export MISS=0.875
export QUAL=50
export MIN_DP=20
export MAX_DP=60
# QUAL/Total DP > 20
# SP > 0.10
# GQ > 40

for NAME in "${VCFLIST[@]}"
do
    echo doing "$NAME" variants
    vcftools --gzvcf "$INDIR"/"$NAME"_snps.vcf.gz --maf $MAF --minQ $QUAL --min-meanDP $MIN_DP --max-meanDP $MAX_DP --minDP $MIN_DP --maxDP $MAX_DP --recode --stdout | bcftools view -i 'QUAL/SMPL_SUM(FORMAT/AD)>20 & FORMAT/SP > 0.10 & FORMAT/GQ > 40' --threads 12 -O v - | vcftools --vcf - --max-missing $MISS --recode --stdout | bcftools sort -m 9500 -O z - > "$NAME"_snps_fil.vcf.gz
done
```

Filtered invariant calls based on stringent parameters.

```bash
# Done via UFRC queue system; see filter_invars.job for more details.
# Resources used: 800 Mb, 2 hrs
module load vcftools/0.1.16
module load bcftools/1.15

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/raw_split"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

export MISS=0.90
export QUAL=50
export MIN_DP=20
export MAX_DP=60
# QUAL/Total DP > 20

for NAME in "${VCFLIST[@]}"
do
    echo doing "$NAME" invariants
    vcftools --gzvcf "$INDIR"/"$NAME"_invar.vcf.gz --minQ $QUAL --min-meanDP $MIN_DP --max-meanDP $MAX_DP --minDP $MIN_DP --maxDP $MAX_DP --remove-indels --recode --stdout | bcftools view -i 'QUAL/SMPL_SUM(FORMAT/AD)>20' --threads 16 -O v - | vcftools --vcf - --max-missing $MISS --recode --stdout | bcftools sort -m 19500 -O z - > "$NAME"_invar_fil.vcf.gz
done
```

Get entry count for each filtered set:

```bash
# Done via UFRC queue system; see fil_var_stats.job for more details.
# Resources used: 2 Mb, 20 sec

module load bcftools/1.15
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_var"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_stats/variants"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    bcftools stats --fasta-ref "$REF_FILE" "$VCF_DIR"/"$NAME"_snps_fil.vcf.gz > "$OUT_DIR"/"$NAME"_snps_fil_stats.txt
done
```

```bash
# Done via UFRC queue system; see fil_invar_stats.job for more details.
# Resources used: 2 Mb, 30 sec

module load bcftools/1.15
REF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_invar"
OUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_stats/invariants"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    bcftools stats --fasta-ref "$REF_FILE" "$VCF_DIR"/"$NAME"_invar_fil.vcf.gz > "$OUT_DIR"/"$NAME"_invar_fil_stats.txt
done
```

## Merge variant and invariant filtered sets
```bash
# Done via UFRC queue system; see merge_fil_vcfs.job for more details.
# Resources used: 5 Gb, 7 min

module load bcftools/1.15

VAR_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_var"
INVAR_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/filtered_invar"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    tabix "$VAR_DIR"/"$NAME"_var_fil.vcf.gz
    tabix "$INVAR_DIR"/"$NAME"_invar_fil.vcf.gz
    bcftools concat --allow-overlaps -O z --threads 12 "$VAR_DIR"/"$NAME"_var_fil.vcf.gz "$INVAR_DIR"/"$NAME"_invar_fil.vcf.gz | bcftools sort -m 4500 -O z - > "$NAME"_fil.vcf.gz
done
```

Moved files to `/orange` drive long-term storage and created symlinks to `/blue` working storage. Merged chromosome VCFs (Do not use this for most analyses, too large!)

```bash
# Done via UFRC queue system; see merge_fil_chr.job for more details.
# Resources used: 2 Mb, 40 sec

module load bcftools/1.15

WDIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/06.snp_filtering"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    tabix "$WDIR"/"$NAME"_fil.vcf.gz
done

ls "$WDIR"/*_fil.vcf.gz > "$WDIR"/vcf_list.txt
bcftools concat -O z --threads 12 -f "$WDIR"/vcf_list.txt -o "$WDIR"/all_fil.vcf.gz
```

Variant/Invariant Counts:

| Chromosome | Num Variants | Num Invariants | **TOTAL**   |
| ---------- | ------------ | -------------- | ----------- |
| Chr01      | 165,262      | 1,259,181      | 1,433,859   |
| Chr02      | 167,481      | 1,368,061      | 1,544,773   |
| Chr03      | 178,061      | 1,180,563      | 1,370,079   |
| Chr04      | 130,871      | 1,002,216      | 1,140,499   |
| Chr05      | 160,414      |   981,011      | 1,151,868   |
| Chr06      | 196,142      | 1,748,766      | 1,955,226   |
| Chr07      | 150,019      | 1,058,566      | 1,217,894   |
| Chr08      | 236,520      | 1,654,303      | 1,904,945   |
| Chr09      | 144,600      | 1,104,510      | 1,257,095   |
| Chr10      | 153,849      | 1,258,706      | 1,420,981   |
| Chr11      | 149,538      | 1,262,938      | 1,420,617   |
| ChrUn      | 7,082        | 34,194         | 41,933      |
| **TOTAL**  | 1,839,839    | 13,913,015     | 15,859,769  |

Filtered SNP set to biallelic SNPs (needed for several analyses) using `vcftools`.

```bash
# Run on UFRC queue system; see get_biallelic.job for more details.
# Resources used: 15 Mb, 7 min

module load vcftools/0.1.16
VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"

vcftools --gzvcf "$VCFDIR"/all_fil.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout | bgzip -c > "$VCFDIR"/all_fil_biallelic.vcf.gz # unzipped this file later for another analysis, so name will be different in files.
```