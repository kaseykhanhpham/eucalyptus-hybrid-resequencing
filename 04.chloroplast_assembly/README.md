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


```