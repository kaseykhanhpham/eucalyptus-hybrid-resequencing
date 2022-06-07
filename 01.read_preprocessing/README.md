# Read Pre-Processing
Housekeeping and quality checks before actually starting any data manipulation

### Read storage and organization
Reads were downloaded from RAPiD and stored in UFRC's `/orange` drive, the longer-term storage drive for data. Most work was done on UFRC's `/blue` drive, for active analysis. A lot of housekeeping throughout will be organizing between the two locations.

**Create symlinks to reads in `\blue` with informative names:**
```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

# Uses information from sample_sequencing_metadata.csv to make informative names
# File addresses are hardcoded into the script
python "$SCRIPT_DIR"/reorg_raw_reads.py
```

### Quality check of raw reads

**Run [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for all files:**
```bash
# Run via job on UFRC, see fastqc_raw.job for details
# Resources used: 357 Mb, 19 hrs

module load fastqc
READ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/raw_reads/"
OUTPUT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/fastqc_raw"

ls "$READ_DIR" | while read NAME; do fastqc -o "$OUTPUT_DIR" "$READ_DIR"/$NAME; done
```

**Summarize FastQC results with [MultiQC](https://multiqc.info/):**
```bash
# Run via job on UFRC, see multiqc_raw.job for details
# Resources used: 5 Gb, 1 min

module load multiqc
FASTQC_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/fastqc_raw"

multiqc "$FASTQC_DIR"
```

### Read trimming and filtering

**Trim sequences with [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic):**
```bash
# Run via job on UFRC, see trimmomatic.job for more details
# Resources: 1 Gb, 4 hrs

module load trimmomatic/0.39
RAW_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/raw_reads"
TRIMMED_DIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/trimmed_reads/"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/"

while read NAME; do trimmomatic PE -threads 11 "$RAW_DIR"/"$NAME"_R1_raw.fastq.gz "$RAW_DIR"/"$NAME"_R2_raw.fastq.gz "$TRIMMED_DIR"/"$NAME"_R1_paired.fq "$TRIMMED_DIR"/"$NAME"_R1_unpaired.fq "$TRIMMED_DIR"/"$NAME"_R2_paired.fq "$TRIMMED_DIR"/"$NAME"_R2_unpaired.fq ILLUMINACLIP:/blue/soltis/kasey.pham/euc_hyb_reseq/reads/euc_adapters.fasta:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36; done < "$LIST_DIR"/seq_ids.txt
```

**Reorganize trimmed reads:**
```bash
INDIR="/orange/soltis/kasey.pham/eucalyptus_hyb_reseq/02.trimmed_reads"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Merge unpaired reads
while read NAME;
do 
    cat "$INDIR"/"$NAME"_R1_unpaired.fq > "$INDIR"/"$NAME"_unpaired.fq
    cat "$INDIR"/"$NAME"_R2_unpaired.fq >> "$INDIR"/"$NAME"_unpaired.fq
done < "$LIST_DIR"/seq_ids.txt
rm *_R*_unpaired.fq

# Make symlinks from /orange to /blue
while read NAME
do
    ln -s "$INDIR"/"$NAME"_R1_paired.fq "$OUTDIR"/"$NAME"_R1_paired_trimmed.fq
    ln -s "$INDIR"/"$NAME"_R2_paired.fq "$OUTDIR"/"$NAME"_R2_paired_trimmed.fq
    ln -s "$INDIR"/"$NAME"_unpaired.fq "$OUTDIR"/"$NAME"_unpaired_trimmed.fq
done < "$LIST_DIR"/seq_ids.txt
```

### Quality check of trimmed reads

**Run `FastQC` for all files:**
```bash
# Run via job on UFRC, see fastqc_trimmed.job for more details
# Resources: 365 Mb, 16 hrs

module load fastqc
READ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
FASTQC_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/fastqc_trimmed"

ls "$READ_DIR"/*.fq | while read NAME; do fastqc -o "$FASTQC_DIR" "$READ_DIR"/$NAME; done
```

**Summarize `FastQC` results with `MultiQC`:**
```bash
# Run via job on UFRC, see multiqc_trimmed.job for details
# Resources used: 5 Gb, 1 min

module load multiqc
FASTQC_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/fastqc_trimmed"

multiqc "$FASTQC_DIR"
```

### Check k-mer distributes of samples

**Run k-mer counting program [`Jellyfish`](https://genome.umd.edu/jellyfish.html) and visualize in [R](https://www.r-project.org/):**

Followed protocol of [this tutorial](https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html) from the Nishiyama lab.
```bash
# Run via job on UFRC, see jellyfish.job for details
# Resources used:  90 Gb, 15 hrs

module load jellyfish/2.3.0
module load R/4.1
READS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

while read NAME
do
    # use a kmer size of 24, a hash size of 5Gb, and 10 threads
    jellyfish count -o "$NAME"_24mer -m 24 -s 5G -t 10 -C "$READS_DIR"/"$NAME"*trimmed.fq
    jellyfish histo -o "$NAME"_24mer.histo "$NAME"_24mer
    # plot k-mer distribution starting at index = 3
    Rscript "$SCRIPT_DIR"/plot_kmer_histo.r "$NAME"_24mer.histo "$NAME"_24mer.png 3
    rm "$NAME"_24mer
done < "$LIST_DIR"/seq_ids.txt
```