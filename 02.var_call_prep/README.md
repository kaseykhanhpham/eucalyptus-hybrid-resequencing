# Preparation for Variant Calling
Steps taken here based loosely on protocol taught in Shin-Han Shiu's Spring 2018 Bioinformatics class at Michigan State University.

### Map to _E. grandis_ reference genome

**Create [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) reference library:**
```bash
module load bowtie2/2.4.2
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"

cd "$REF_DIR"
# Create a bowtie2 index with the name ASM1654582v1
bowtie2-build GCF_016545825.1_ASM1654582v1_genomic.fna ASM1654582v1
```

**Map reads in bowtie2:**
```bash
# Run via job on UFRC, see map_reads.job for details
# Resources used: 953 Mb, 40 hrs

module load bowtie2/2.4.2
READS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/map_reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

export BOWTIE2_INDEXES="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"


while read NAME
do 
    bowtie2 -p 16 -q --end-to-end -x ASM1654582v1 -1 "$READS_DIR"/"$NAME"_R1_paired_trimmed.fq -2 "$READS_DIR"/"$NAME"_R2_paired_trimmed.fq -S "$OUTDIR"/"$NAME"_mapped.sam --no-unal
done < "$LIST_DIR"/seq_ids.txt
```

**Convert files to BAM format using [samtools](https://github.com/samtools/samtools):**
```bash
# Run via job on UFRC, see sam2bam.job for details
# Resources used: 27 Mb, 2 hrs

module load samtools/1.12
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do 
    samtools view -b -@ 8 -o "$NAME"_mapped.bam "$NAME"_mapped.sam
done < "$LIST_DIR"/seq_ids.txt
```

### Add readgroups
Assigning readgroups is necessary because variant callers use information from what lane and run a sample was on to correct for calling bias. I based the information included on the [GATK recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) and [Sentieon recommendations](https://support.sentieon.com/appnotes/read_groups/).

**Make job to assign the right readgroups to each BAM file**

I manually made a table (readgroup_table.txt) with the information for readgroup assignment from my [sample sequencing metadata](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/00.metadata/03.seq_analysis/sample_sequencing_metadata_all.xlsx). The python script I wrote to make the job has a specific file format hard-coded in, but it can be tweaked to fit other formats.

```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/map_reads"
IN_SUFFIX="_mapped.bam"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads"
OUT_SUFFIX="_wRG.bam"
UFRC_ACCT="soltis"
UFRC_EMAIL="kasey.pham@ufl.edu"

python "$SCRIPT_DIR"/make_rg_job.py readgroup_table.txt "$INDIR" "$IN_SUFFIX" "$OUTDIR" "$OUT_SUFFIX" "$UFRC_ACCT" "$UFRC_EMAIL"
```

**Add readgroups using [Picard](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-):**

```bash
# Run via job on UFRC, see add_readgroups.job for details
# Resources used: 10 Gb, 12 hrs
# Commands for each file unique, only showing one example below. Refer to linked resources above for explanation of what each flag means.

module load picard/2.25.5
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/map_reads"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/01.readgroups"

picard AddOrReplaceReadGroups -I "$INDIR"/S320_mapped.bam -O "$OUTDIR"/S320_wRG.bam -LB WA01.leaf -PL ILLUMINA -PU H3K2CDSX2.4.97 -SM WA01 -ID H3K2CDSX2.4.97
```

### Fix mate pairs using samtools

```bash
# Run via job on UFRC, see fixmate.job for details
# Resources used: 15 Mb, 2 hrs

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/01.readgroups"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/02.fixmate"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    echo doing "$NAME" 1>&2 
    samtools fixmate -m -O bam -@ 7 "$INDIR"/"$NAME"_wRG.bam "$OUTDIR"/"$NAME"_fixed.bam
done < "$LIST_DIR"/seq_ids.txt
```

### Sort reads using samtools

```bash
# Run via job on UFRC, see sortbam.job for details
# Resources used: 8 Gb, 4 hrs

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/02.fixmate"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/03.sort"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Sort using a hash memory of 1Gb and 7 threads
while read NAME
do 
    samtools sort -m 1G -o "$OUTDIR"/"$NAME"_sorted.bam -O bam -@ 7 "$INDIR"/"$NAME"_fixed.bam
done < "$LIST_DIR"/seq_ids.txt
```

### Mark duplicates

**Mark duplicate reads in samtools:**

```bash
# Run via job on UFRC, see markdupes.job for details
# Resources used: 158 Mb, 2 hrs

module load samtools/1.12
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/03.sort"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/04.markdup"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Mark duplicates with an optical distance of 2500
while read NAME
do 
    samtools markdup -d 2500 -@ 7 "$INDIR"/"$NAME"_sorted.bam "$OUTDIR"/"$NAME"_marked.bam
done < "$LIST_DIR"/seq_ids.txt
```

**Index BAM files:**

```bash
# Run via job on UFRC, see index.job for details
# Resources used: 10 Mb, 13 min

module load samtools/1.12
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

module load samtools/1.12

while read NAME
do 
    samtools index -b -@ 7 "$NAME"_marked.bam "$NAME"_marked.bai
done < "$LIST_DIR"/seq_ids.txt
```

## Calculate coverage

Wrote `R` script to generate job that calculates depth for each sample at `make_depth_job.r`. Values for directories, etc. are currently hardcoded into the script.

```bash
module load R

Rscript make_depth_job.r
```

**Calculate depth distribution for each sample using samtools:**
```bash
# Run via job on UFRC, see index.job for details
# Resources used: 8 Mb, 5 hrs
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
# Resources used: 8 Mb, 5 hrs

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

### Housekeeping; delete intermediate files:
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
