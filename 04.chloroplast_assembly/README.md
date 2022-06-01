# Chloroplast Assembly

### Run `FastPlast`

**Create job to run `FastPlast`:**
Wrote python script to generate a job file which runs [`FastPlast`](https://github.com/mrmckain/Fast-Plast) on all samples, using trimmed reads from both sequencing runs for each sample:

```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPT_DIR"/generate_fp_job.py
```

**Run `FastPlast`:**
```bash
# Run via job on UFRC, see fastplast.job for details
# Showing command for 1 out of 40 samples as commands are essentially the same for all, save for sample IDs.
# Resources used: 

READ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/organelle"

perl $HPC_FASTPLAST_DIR/Fast-Plast/fast-plast.pl -1 "$READ_DIR"/S320_R1_paired.fq,"$READ_DIR"/S1_R1_paired.fq -2 "$READ_DIR"/S320_R2_paired.fq,"$READ_DIR"/S1_R2_paired.fq -n WA01 --subsample 30000000 --threads 7 --user_bowtie "REF_DIR"/AY780259.1 --clean deep --skip trim --coverage_analysis --min_coverage 10
```
