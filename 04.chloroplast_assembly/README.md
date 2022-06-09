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
# Run via job on UFRC, see fastplast1.job, fastplast2.job, fastplast3.job, fastplast4.job for details
# Showing command for 1 out of 40 samples as commands are essentially the same for all, save for sample IDs.
# Resources used: 5 Gb, 9 hrs
# Bowtie2 always throws an error and coredumps at some point in mapping the subset of reads, but I am getting good looking chloroplasts anyways and it is not because my trimmed reads are corrupted (they map with Bowtie2 to the reference plastome just fine)... So I am not going to try to fix it.

READ_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/reads/trimmed_reads"
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/organelle"

perl $HPC_FASTPLAST_DIR/Fast-Plast/fast-plast.pl -1 "$READ_DIR"/S320_R1_paired.fq,"$READ_DIR"/S1_R1_paired.fq -2 "$READ_DIR"/S320_R2_paired.fq,"$READ_DIR"/S1_R2_paired.fq -n WA01 --subsample 30000000 --threads 7 --user_bowtie "REF_DIR"/AY780259.1 --clean deep --skip trim --coverage_analysis --min_coverage 10
```

**Housekeeping:**
```bash
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WA01 WA02 WA03 WA04 WA05 WB01 WB02 WB03 WB04 WB05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done

SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
for NAME in WC01 WC02 WC03 WC04 WC05 WD01 WD02 WD03 WD04 WD05
do
    "$SCRIPTS_DIR"/clean_fastplast_dir.sh $NAME
done
```