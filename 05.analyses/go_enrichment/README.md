# Gene Ontology Enrichment of Genome Outlier Windows

## Get genes within pi/dxy outlier windows
Remove headers from Dxy/Pi overlap files to make them BED format

```bash
INDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/no_outgroup/5k_windows/pi_samplewise/common_windows/Dxy_Pi"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

while read NAME;
do
    tail -n +2 "$INDIR"/common_5k_Dxy.tab_5k_Pi_"$NAME".tab| awk '{gsub(" ","\t"); print}' > common_Dxy_Pi_"$NAME".bed
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Get genes which fall within outlier windows using `BEDTools`.

```bash
# Performed in UFRC queue system. See bedtools.job for more details.
# Resources:

module load bedtools/2.30.0

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/no_outgroup/go_enrichment"
WINDOWS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/no_outgroup/go_enrichment/dxy_pi_windows"
ANNOT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

while read NAME
do
    bedtools intersect -a "$WINDOWS_DIR"/common_Dxy_Pi_"$NAME".bed -b "$ANNOT_DIR"/EGLOB-X46.v1.0.annotation.gff -wb > common_Dxy_Pi_genes_"$NAME".tab
done < "$LIST_DIR"/Eglobulus_MR.txt
```

Process output of `BEDTools` to write lists of genes within the outlier windows.

```bash
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
module load python/3.8

while read NAME
do
    python parse_intersects.py ./dxy_pi_windows/common_Dxy_Pi_"$NAME".bed common_Dxy_Pi_genes_"$NAME".tab
done < "$LIST_DIR"/Eglobulus_MR.txt
```

## Get gene ontology terms for _E. globulus_ genome annotation

Preliminary annotations for the _E. globulus_ draft genome were transferred over from the _E. grandis_ reference genome using `Liftoff` by Jakob Butler. Used [`GOMAP`](https://dill-picl.org/projects/gomap/) to assign gene ontology terms to gene models.

First extracted gene models from the genome sequence using `BEDTools`.

```bash
module load emboss/6.6.0
module load bedtools/2.30.0
GENOME_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

# Get CDS entries from bedtools and indicate which blocks go together
Rscript "$SCRIPT_DIR"/get_annot_bed.r "$GENOME_DIR"/EGLOB-X46.v1.0.annotation.gff EGLOB-X46.v1.0.annotation_cds.bed12

# Extract concatenated CDS sequences from genome using annotation
bedtools getfasta -fi "$GENOME_DIR"/EGLOB-X46.v1.0.fa -bed EGLOB-X46.v1.0.annotation_cds.bed12 -split -nameOnly | fold -w 80 > EGLOB-X46.v1.0.genes.fa

# translate to AA sequence
transeq EGLOB-X46.v1.0.genes.fa EGLOB-X46.v1.0.aa.fa


```