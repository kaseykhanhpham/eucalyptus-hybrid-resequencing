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