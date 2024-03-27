# GO Enrichment Analysis
Performed Gene Ontology enrichment analysis on regions of interest (ROIs) identified in genome_scan/chr_plots step. Used `LiftOff` annotations of _E. globulus_ X46 genome; these are ported from the _E. grandis_ annotation, so they should be regarded as a very rough draft and taken with a grain of salt.

First, converted _E. globulus_ X46 annotation to BED format in `R v.4.2`.
```R
# REDUNDANT?
annot_name <- "/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.annotation.gff" # could do it with my re-run of LiftOff instead if I want the E. grandis gene IDs
annot_tab <- read.table(annot_name, header = FALSE, sep = "\t")
# subset to just genes
annot_tab <- annot_tab[which(annot_tab$V3 == "gene"),]
# extract source gene ID
notes_id <- unlist(lapply(strsplit(annot_tab$V9, ";"), function(x) x[1]))
id <- unlist(sapply(notes_id, function(x) gsub("ID=", "", x, fixed = TRUE)))
# build final output BED
out_tab <- data.frame(chrom = annot_tab$V1, start = annot_tab$V4, end = annot_tab$V5, name = id)
# export
out_name <- "/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.annotation.bed"
write.table(out_tab, out_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
```

## Species differences
Retrieved genes with at least 75% overlap with species difference regions of interest using `BEDTools`.

```bash
module load bedtools
ANNOT_BED="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.annotation_cds.bed12"
FDIFFS_BED="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/chr_plots/fdiff_files/Eglob_all_fdiffs.bed"
GENOME="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.fa"

# get overlapping annotated genes and the number of base pairs overlapping with the species diffs BED file
bedtools intersect -a "$ANNOT_BED" -b "$FDIFFS_BED" -u -f 0.75 -header > Eglob_all_fdiffs_genes.bed
# extract genes from genome FASTA
bedtools getfasta -fi "$GENOME" -bed Eglob_all_fdiffs_genes.bed -nameOnly > Eglob_all_fdiffs_genes.fas
```

Retrieved 312 genes with >75% overlap with species differences (FST between _E. cordata_ and _E. globulus_ in 0.95 percentile).

Ran [`GOMAP Singularity`](https://github.com/Dill-PICL/GOMAP-singularity/tree/master) for species difference genes. See fdiffs_gomapX.job for run details, but the execution is done through wrapper scripts, so there isn't much to it.

## Recombination suppression
Checked how much overlap there was with genes to begin with.

```bash

```