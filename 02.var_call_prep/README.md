# Preparation for Variant Calling
Steps taken here based loosely on protocol taught in Shin-Han Shiu's Spring 2018 Bioinformatics class at Michigan State University.

### Map to _E. globulus_ reference genome

**Create [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) reference library:**
```bash
module load bwa-mem2/2.2.1
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"

cd "$REF_DIR"
# Create a bwa-mem index with the name EGLOB-X46.v1
bwa-mem2 index -p EGLOB-X46.v1 EGLOB-X46.v1.0.fa
```

**Map reads in `bwa-mem2`:**
```bash
# Run via job on UFRC, see map_reads.job for details
# Resources used: 8 Gb, 31 hrs

module load bwa-mem2/2.2.1
READS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/01.map_reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
INDEX_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"

while read NAME
do
    echo doing "$NAME"
    bwa-mem2 mem -t 16 -o "$OUTDIR"/"$NAME"_mapped.sam "$INDEX_DIR"/EGLOB-X46.v1 "$READS_DIR"/"$NAME"_R1_paired.fq "$READS_DIR"/"$NAME"_R2_paired.fq
done < "$LIST_DIR"/seq_ids.txt
```

**Get mapping rate:**
```bash
# Run via job on UFRC, see mapping_stats.job for details
# Resources used: 20 Mb, 16 min
module load samtools/1.12
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    echo doing "$NAME"
    samtools flagstats -@ 8 "$NAME"_mapped.sam
done < "$LIST_DIR"/seq_ids.txt
```

Saved a CSV file with the sequence IDs and their resulting mapping rate log. Checked mapping rate between species on local install of `R`:

```R
# set directories and input files
working_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/02.var_call_prep"
metadata_dir <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis"

mapping_rate_file <- paste(working_dir, "bwa_mapping_rate.csv", sep = "/")
seq_to_sample_file <- paste(metadata_dir, "sample_sequencing_metadata.csv", sep = "/")
sample_to_spp_file <- paste(metadata_dir, "sample_spp_table.csv", sep = "/")

setwd(working_dir)

# import input files
mapping_rate <- read.csv(mapping_rate_file)
seq_to_sample <- read.csv(seq_to_sample_file, as.is = TRUE)
sample_to_spp <- read.csv(sample_to_spp_file, as.is = TRUE)

# format metadata
seq_to_sample <- unique(seq_to_sample[,c("Sample", "RunSample")])

# associate sample ids with species
sample_ids <- seq_to_sample[match(mapping_rate$seq_id, seq_to_sample$RunSample), "Sample"]
spp_ids <- sample_to_spp[match(sample_ids, sample_to_spp$RAPiD_ID), "Taxon"]
mapping_rate <- cbind(mapping_rate, spp_ids = as.factor(spp_ids))

## Testing mapping rate difference between species
# ANOVA test
map_anova <- aov(map_rate ~ spp_ids, data = mapping_rate) # P = 0.165

# Tukey's pairwise comparisons
mean(mapping_rate[which(mapping_rate$spp_ids == "glob_pure"),"map_rate"]) # 0.9856
mean(mapping_rate[which(mapping_rate$spp_ids == "glob_MR"),"map_rate"]) # 0.9864
mean(mapping_rate[which(mapping_rate$spp_ids == "cord_MR"),"map_rate"]) # 0.9838

TukeyHSD(map_anova, conf.level = 0.95) # no significant difference between any population comparisons

## Testing pairing rate difference between species
# ANOVA test
pair_anova <- aov(pair_rate ~ spp_ids, data = mapping_rate) # P = 0.0593

# Tukey's pairwise comparisons
mean(mapping_rate[which(mapping_rate$spp_ids == "glob_pure"),"pair_rate"]) # 0.9856
mean(mapping_rate[which(mapping_rate$spp_ids == "glob_MR"),"pair_rate"]) # 0.9864
mean(mapping_rate[which(mapping_rate$spp_ids == "cord_MR"),"pair_rate"]) # 0.9838

```
### Mapping Rate Stats
ANOVA results:

| Component | df | SS         | MSE        | F     | P-value |
| --------- | -- | ---------- | ---------- | ----- | ------- |
| species   | 2  | 0.0000892  | 0.00004460 | 1.846 | 0.165   |
| residuals | 77 | 0.001861   | 0.00002417 | -     | -       |

Tukey results:
| Comparison                             | Difference | Conf. Int. Lower | Conf. Int. Upper | P-value (adj) |
| -------------------------------------- | ---------- | ---------------- | ---------------- | ------------- |
| introgressed E. globulus to E. cordata | 0.002585   | -0.0006324       | 0.005802         | 0.140         |
| reference E. globulus to E. cordata    | 0.001810   | -0.001905        | 0.005525         | 0.478         |
| introgressed E. globulus to reference E. globulus | -0.003992 | -0.00776 | 0.002442       | 0.833         |

### Pairing Rate Stats
ANOVA results:

| Component | df | SS      | MSE       | F     | P-value |
| --------- | -- | ------- | --------- | ----- | ------- |
| species   | 2  | 0.00324 | 0.001619  | 2.931 | 0.0593  |
| residuals | 77 | 0.04253 | 0.0005524 | -     | -       |

Tukey results:
| Comparison                             | Difference | Conf. Int. Lower | Conf. Int. Upper | P-value (adj) |
| -------------------------------------- | ---------- | ---------------- | ---------------- | ------------- |
| introgressed E. globulus to E. cordata | 0.01467    | -0.0007049       | 0.03006          | 0.0646        |
| reference E. globulus to E. cordata    | 0.004845   | -0.01292         | 0.02261          | 0.792         |
| introgressed E. globulus to reference E. globulus | -0.009832 | -0.02521 | 0.005550       | 0.284         |

None of the differences in mapping or pairing rate was large or significant enough for me to be particularly concerned. Moved on with analysis.

**Convert files to `BAM` format using [`samtools`](https://github.com/samtools/samtools):**
NOTE (TO REMOVE LATER): can re-run this job as-is
```bash
# Run via job on UFRC, see sam2bam.job for details
# Resources used: 31 Mb, 1.5 hrs

module load samtools/1.12
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do 
    samtools view -b -@ 8 -o "$NAME"_mapped.bam "$NAME"_mapped.sam
done < "$LIST_DIR"/seq_ids.txt
```

### Add readgroups
Assigning readgroups is necessary because variant callers use information from what lane and run a sample was on to correct for calling bias. I based the information included on the [`GATK` recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) and [Sentieon recommendations](https://support.sentieon.com/appnotes/read_groups/).

**Make job to assign the right readgroups to each `BAM` file**

I manually made a table (readgroup_table.txt) with the information for readgroup assignment from my [sample sequencing metadata](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/00.metadata/03.seq_analysis/sample_sequencing_metadata_all.xlsx). The python script I wrote to make the job has a specific file format hard-coded in, but it can be tweaked to fit other formats.

```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/01.map_reads"
IN_SUFFIX="_mapped.bam"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/01.readgroups"
OUT_SUFFIX="_wRG.bam"
UFRC_ACCT="soltis"
UFRC_EMAIL="kasey.pham@ufl.edu"

python "$SCRIPT_DIR"/make_rg_job.py readgroup_table.txt "$INDIR" "$IN_SUFFIX" "$OUTDIR" "$OUT_SUFFIX" "$UFRC_ACCT" "$UFRC_EMAIL"
```

**Add readgroups using [`Picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-):**

```bash
# Run via job on UFRC, see add_readgroups.job for details
# Resources used: 10 Gb, 19 hrs
# Commands for each file unique, only showing one example below. Refer to linked resources above for explanation of what each flag means.

module load picard/2.25.5
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/01.map_reads"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/01.readgroups"

picard AddOrReplaceReadGroups -I "$INDIR"/S320_mapped.bam -O "$OUTDIR"/S320_wRG.bam -LB WA01.leaf -PL ILLUMINA -PU H3K2CDSX2.4.97 -SM WA01 -ID H3K2CDSX2.4.97
```

### Fix mate pairs using `samtools`

```bash
# Run via job on UFRC, see fixmate.job for details
# Resources used: 15 Mb, 2 hrs

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/01.readgroups"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/02.fixmate"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    echo doing "$NAME" 1>&2 
    samtools fixmate -m -O bam -@ 7 "$INDIR"/"$NAME"_wRG.bam "$OUTDIR"/"$NAME"_fixed.bam
done < "$LIST_DIR"/seq_ids.txt
```

### Sort reads using `samtools`

```bash
# Run via job on UFRC, see sortbam.job for details
# Resources used: 7.5 Gb, 5 hrs

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/02.fixmate"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/03.sort"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Sort using a hash memory of 1Gb and 7 threads
while read NAME
do 
    samtools sort -m 1G -o "$OUTDIR"/"$NAME"_sorted.bam -O bam -@ 7 "$INDIR"/"$NAME"_fixed.bam
done < "$LIST_DIR"/seq_ids.txt
```

### Mark duplicates

**Mark duplicate reads in `samtools`:**

```bash
# Run via job on UFRC, see markdupes.job for details
# Resources used: 350 Mb, 2 hrs

module load samtools/1.12

INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/03.sort"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
LISTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    samtools markdup -d 2500 -@ 7 "$INDIR"/"$NAME"_sorted.bam "$OUTDIR"/"$NAME"_marked.bam
done < "$LISTDIR"/seq_ids.txt
```

**Index `BAM` files:**

```bash
# Run via job on UFRC, see index.job for details
# Resources used: 9 Mb, 10 min

module load samtools/1.12

LISTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    samtools index -b -@ 7 "$NAME"_marked.bam "$NAME"_marked.bai
done < "$LISTDIR"/seq_ids.txt
```

## Calculate coverage

Wrote `R` script to generate job that calculates depth for each sample at `make_depth_job.r`. Values for directories, etc. are currently hardcoded into the script.

```bash
module load R/4.2
SCRIPTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

Rscript "$SCRIPTDIR"/make_depth_job.r
```

**Calculate depth distribution for each sample using `samtools`:**
```bash
# Run via job on UFRC, see index.job for details
# Resources used: 26 Mb, 6 hrs
# Commands for each file unique because runfiles for each sample must be matched, only showing one example below.

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/04.markdup"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/05.depth"

samtools depth -o "$OUTDIR"/WA01_cover.txt "$INDIR"/S320_marked.bam "$INDIR"/S1_marked.bam
```

**Plot coverage and calculate average:**

Used custom python script for this; it calculates the average fine, but the plotting looks a bit funny. I never quite worked out what the issue was.
```bash
# Run via job on UFRC, see coverage.job for details
# Resources used: 10 hrs, 70 Mb

module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# create blank file in which to dump average coverage calculations
touch average_coverage.txt

# calculate average from samtools output and plot average coverage across all chromosomes in sliding window
while read NAME
do 
    python "$SCRIPT_DIR"/plot_coverage.py "$NAME"_cover.txt "$NAME" average_coverage.txt
done < "$LIST_DIR"/sample_ids.txt

```

## Housekeeping

**Delete intermediate files:**
```bash
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/map_reads"
RG_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/01.readgroups"
FIXMATE_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/02.fixmate"
SORT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/03.sort"
DEPTH_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/05.depth"

rm "$MAP_DIR"/*.sam
rm "$RG_DIR"/*_wRG.bam
rm "$FIXMATE_DIR"/*_fixed.bam
rm "$SORT_DIR"/*_sorted.bam
rm "$DEPTH_DIR"/*_cover.txt
```

**Reorganize final files for SNP calling:**
```bash
READ_INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/02.process_reads/04.markdup"
READ_OUTDIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/04.processed_reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Move to file storage
mv "$READ_INDIR"/*.bam "$READ_OUTDIR"

# Create symlinks to working directory
while read NAME; do ln -s "$READ_OUTDIR"/"$NAME"_marked.bam "$READ_INDIR"/"$NAME"_marked.bam; done < "$LIST_DIR"/seq_id.txt
```