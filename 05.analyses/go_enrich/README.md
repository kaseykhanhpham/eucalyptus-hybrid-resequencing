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
bedtools getfasta -fi "$GENOME" -bed Eglob_all_fdiffs_genes.bed -nameOnly > recip_blast/Eglob_all_fdiffs_genes.fas
```

Retrieved 312 genes with >75% overlap with species differences (FST between _E. cordata_ and _E. globulus_ in 0.95 percentile).

Reciprocal `BLAST` between _E. grandis_ genome annotation and species difference genes.

```bash
module load ncbi_blast/2.15.0
cd /blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/annotation
makeblastdb -in Egrandis_297_v2.0.cds.fa -parse_seqids -dbtype nucl
cd /blue/soltis/kasey.pham/euc_hyb_reseq/analyses/go_enrich/fdiffs/recip_blast
makeblastdb -in Eglob_all_fdiffs_genes.fas -parse_seqids -dbtype nucl
```

```bash
# Run on UFRC queue system; see fdiffs_blast.job for more details.
# Resources used: 15 sec, 1 Gb

module load ncbi_blast/2.15.0
FDIFFS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/go_enrich/fdiffs/recip_blast"
REFS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/annotation"
export BLAST_DB="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/annotation:/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/go_enrich/fdiffs/recip_blast"

blastn -query "$FDIFFS_DIR"/Eglob_all_fdiffs_genes.fas -db Egrandis_297_v2.0.cds.fa -outfmt 6 -evalue 0.01 -num_threads 8 -out Eglob_fdiffs_to_AUSX01_cds_BLAST.tab

blastn -query "$REFS_DIR"/Egrandis_297_v2.0.cds.fa -db Eglob_all_fdiffs_genes.fas -outfmt 6 -evalue 0.01 -num_threads 8 -out AUSX01_cds_to_Eglob_fdiffs_BLAST.tab
```

Retrieved reciprocal BLAST hits from the two searches. Extracted corresponding homologous protein sequence in _E. grandis_ reference genome annotation using `R` script.

```bash
module load R/4.2
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/annotation"

# Get homologous E. grandis transcripts from reciprocal best matches
Rscript "$SCRIPT_DIR"/get_rblast_matches.r -i Eglob_fdiffs_to_AUSX01_cds_BLAST.tab -j AUSX01_cds_to_Eglob_fdiffs_BLAST.tab -l ../Eglob_fdiffs_AUSX01_homol_prot_list.txt -t ../Eglob_fdiffs_AUSX01_homol_tab.txt -n 2

# Get E. grandis protein seqs
cd ..
Rscript "$SCRIPT_DIR"/extract_fas.r -i "$REF_DIR"/Egrandis_297_v2.0.protein.fa -l Eglob_fdiffs_AUSX01_homol_prot_list.txt -o Eglob_fdiffs_AUSX01_homol_prot.fas
```

Ran `GO-MAP` on _E. grandis_ homologous protein sequences for approximate GO enrichment.

## Recombination
Checked how much overlap there was with genes to begin with.

```bash

```

## Introgression

```bash
```

## Selection

```bash
```