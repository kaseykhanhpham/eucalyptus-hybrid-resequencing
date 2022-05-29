# Variant Calling

### Housekeeping:
```bash
READ_INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/process_reads/04.markdup"
READ_OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes/reads"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"
VC_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes"

# create symlinks to processed reads in variant calling directory
while read NAME; do ln -s "$READ_INDIR"/"$NAME"_marked.bam "$READ_OUTDIR"/"$NAME"_marked.bam; done < "$LIST_DIR"/sample_ids.txt
# create list of BAM inputs for variant caller
ls -d "$READ_OUTDIR"/* > "$VC_DIR"/bam_inputs.txt
```

### Run Variant Caller
I am using the variant caller [FreeBayes](https://github.com/freebayes/freebayes), which estimates internal population parameters for the individuals it is provided with using Bayesian inference. This is probably a good choice for my dataset given that I don't have a set of validated SNPs for training a variant caller like GATK. I followed recommendations on FreeBayes' Github page for parameter selection.

I ran FreeBayes on each chromosome separately, with the same parameters each time. See the job files named `fb_NC0526**.job` for each set of individual parameters.

**Command for chromosome 1:**

```bash
# Run via job on UFRC, see fb_NC052612.job for details
# Resources used: 3 Gb, 7 days 
# (range: 2.5 Gb - 4.5 Gb, 2 days - 11 days)

module load freebayes/1.3.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes/max_dp_3000"

# Run FreeBayes on all samples at once against E. grandis reference genome, call sites with a maximum depth of 3000
freebayes -L "$LIST_DIR"/bam_inputs.txt -f "$REF_DIR"/GCF_016545825.1_ASM1654582v1_genomic.fna -r NC_052612.1 -v "$OUTDIR"/NC_052612.vcf -g 3000
```

**Command for unanchored contigs:**

All unanchored contigs in reference genome were called together. I manually made a BED file of the names of each contig as input for FreeBayes.

```bash
# Run via job on UFRC, see fb_contigs.job for details
# Resources used: 4.1 Gb, 2 days 

module load freebayes/1.3.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/ASM1654582v1"
OUTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes/max_dp_3000"


freebayes -L "$LISTDIR"/bam_inputs.txt -f "$REF_DIR"/GCF_016545825.1_ASM1654582v1_genomic.fna -t "$LIST_DIR"/ASM1654582_unanchored_contigs.bed  -v "$OUTDIR"/unanchored_contigs.vcf
```

**Merge VCF files for each chromosome using [picard](https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard):**

```bash
# Run via job on UFRC, see merge_vcfs.job for details
# Resources used: 10.8 Gb, 2 hrs 

module load picard/2.25.5
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes/max_dp_3000"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/freebayes"

ls "$VCF_DIR"/*.vcf > "$LIST_DIR"/vcfs.list
picard MergeVcfs -I "$LIST_DIR"/vcfs.list -O all_to_ASM1654582v1.vcf
```

