# Sliding Window Genome Statistics

## Nucleotide Diversity and Absolute Divergence

Used [`pixy`](https://pixy.readthedocs.io/en/latest) to calculate pi (nucleotide diversity) and dxy (absolute divergence) across the genome.

```bash
# Done on UFRC queue system; see pixy.job for more details.
# Resources used:

module load conda 

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy"

conda activate "$ENV_DIR"/pixy

pixy --stats pi fst dxy --vcf "$WDIR"/meehan_all_fil_maf0.00_snps_noout.vcf --populations "$WDIR"/populations.txt --window_size 5000 --n_cores 12 --output_prefix meehan_maf0.00 --chromosomes Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11 --fst_type wc
```