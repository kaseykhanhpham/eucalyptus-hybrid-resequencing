# Annotation of X46 Genome Using LiftOff

LiftOff was run with polishing to transfer annotation from the _E. grandis_ AUSX01v2 genome annotation to the _E. globulus_ X46 assembly. This had been done before; re-ran with polishing to correct exon boundaries.

```bash
# Run on UFRC queue system; see liftoff.job for more details.
# Resources used: 18 Gb, 5 min

module load liftoff/1.4.2
TARGET_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
ASSEMB_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/assembly"
ANNOT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/AUSX01_v2/annotation"

liftoff -g "$ANNOT_DIR"/Egrandis_297_v2.0.gene.gff3 -o Eglob_X46_liftoff_20230731.gff3 -u Eglob_X46_liftoff_20230731_unmapped.txt -exclude_partial -a 0.7 -s 0.7 -flank 0.1 -chroms chr_X46_to_AUSX01.txt -unplaced Eglob_X46_liftoff_20230731_unplaced.txt -copies -sc 0.975 -p 9 -polish "$TARGET_DIR"/EGLOB-X46.v1.0.fa "$ASSEMB_DIR"/Egrandis_297_v2.0.fa
```
Separated the resulting GFF3 file into just gene and just CDS GFF files in `python 3`.

```python
infile_name = "Eglob_X46_liftoff_20230731.gff3_polished"
gene_name = "Eglob_X46_liftoff_20230731_genes.gff3_polished"
cds_name = "Eglob_X46_liftoff_20230731_cds.gff3_polished"

infile = open(infile_name, "r")
out_genes = open(gene_name, "w")
out_cds = open(cds_name, "w")

for line in infile:
    if line.startswith("#"):
        out_genes.write(line)
        out_cds.write(line)
    elif line.strip().split()[2] == "gene":
        out_genes.write(line)
    elif line.strip().split()[2] == "CDS":
        out_cds.write(line)

infile.close()
out_genes.close()
out_cds.close()
```

Extracted genes, CDS, and translated proteins from the LiftOff annotation using `BEDTools`. 
```bash
module load bedtools/2.30.0
SEQ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"

bedtools getfasta -fo Eglob_X46_liftoff_20230731_genes.fa -fi "$SEQ_DIR"/EGLOB-X46.v1.0.fa -bed Eglob_X46_liftoff_20230731_genes.gff3_polished
bedtools getfasta -fo Eglob_X46_liftoff_20230731_cds.fa -fi "$SEQ_DIR"/EGLOB-X46.v1.0.fa -bed Eglob_X46_liftoff_20230731_cds.gff3_polished
```

Translated CDS to protein sequences. (Still need to concatenate CDS from same gene...)
```bash
module load emboss/6.6.0

transeq Eglob_X46_liftoff_20230731_cds.fa Eglob_X46_liftoff_20230731_pep.fa
```