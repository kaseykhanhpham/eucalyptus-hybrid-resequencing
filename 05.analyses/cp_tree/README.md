# Chloroplast Haplotype Analysis

## Get chloroplast alignment

**Create FASTA files with all chloroplast region assemblies and reformat:**
Wrote python scripts to concatenate chloroplast assemblies from each sample into separate FASTA files for each region (files currently hardcoded!) and reformat FASTA files to split sequences into multiple lines.

```bash
module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPTS_DIR"/make_cp_region_fas.py
# split each sequence into multiple lines
python "$SCRIPTS_DIR"/reformat_fasta.py 80 irb_temp.fas irb_assemblies.fas
python "$SCRIPTS_DIR"/reformat_fasta.py 80 lsc_temp.fas lsc_assemblies.fas
python "$SCRIPTS_DIR"/reformat_fasta.py 80 ssc_temp.fas ssc_assemblies.fas
rm *_temp.fas
```

**Align chloroplast regions:**
```bash
# Run via job queue on UFRC, see mafft_ssc.job, mafft_irb.job, mafft_lsc.job for more details.
# Resources used: ~500 Mb, 11 min maximum.
# One example of the three regions shown below.

module load mafft/7.490

mafft --auto --thread 16 lsc_assemblies.fas > lsc_assemblies_aligned.fas
```

Sequences looked good on manual inspection, decided not to mess with them and move on to phylogenetic inference.

**Concatenate regions:**
```bash
module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPTS_DIR"/concatenate_fasta.py irb_assemblies_aligned.fas,lsc_assemblies_aligned.fas,ssc_assemblies_aligned.fas concatenated_cp_aligned.fas
```

## Infer chloroplast tree
