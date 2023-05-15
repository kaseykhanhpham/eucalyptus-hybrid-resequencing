# HMM Inference of Local Ancestry
Thanks to Shelley Sianta for suggesting this analysis.

## Prepare inputs
### Genotype count input file

First retrieved allele counts for each sample using `VCFTools`.

```bash
# Performed in UF HPC queue system; see gt_counts.job for more details.
# Resources used: 4 Mb, 40 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --vcf "$INFILE" --out cord_counts --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --counts

# E. globulus ref
vcftools --vcf "$INFILE" --out glob_counts --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --counts

# Introgressants
while read NAME
do
    vcftools --vcf "$INFILE" --out "$NAME"_counts --indv "$NAME" --min-alleles 1 --max-alleles 2 --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt

# Removed unanchored contigs from counts files because they are not in genetic map
ls *.count | while read NAME; do head -n 9703613 "$NAME" > "$NAME".tmp; done
remove *.count
rename count.tmp count *
```

Used custom `R` and `python` scripts to generate input for `Ancestry_HMM` from `VCFTools` counts.
```bash
# Performed in UF HPC queue system; see get_ahmm_in.job for more details
# Resources used:

module load R/4.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_counts.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt coverage.txt sample_files.txt gen_dists.tab meehan_all_ahmm_in.tab
```

