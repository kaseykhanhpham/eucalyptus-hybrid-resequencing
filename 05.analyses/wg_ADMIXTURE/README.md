# Whole-genome `ADMIXTURE`
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/ADMIXTURE/).
Uses the maximum-likelihood population structure program [`ADMIXTURE`](http://dalexander.github.io/admixture/).

## `ADMIXTURE` on all samples

### Create `BIM` file with variants for `ADMIXTURE`
This analysis needs a set of unlinked SNPs as input. I already generated a set of variants pruned for linkage for the PCA analysis, see [that section](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/pca) for more details.

```bash
module load plink/1.90b3.39
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05"
PRUNED_VARS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca"

plink --vcf "$VCF_DIR"/all_to_ASM1654582_fil_maf0.05.vcf --extract "$PRUNED_VARS_DIR"/all_maf0.05.prune.in --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out all_to_ASM1654582_fil_maf0.05
awk '{$1="0";print $0}' all_to_ASM1654582_fil_maf0.05.bim > all_to_ASM1654582_fil_maf0.05.bim.tmp
mv all_to_ASM1654582_fil_maf0.05.bim.tmp all_to_ASM1654582_fil_maf0.05.bim
```

### Run `ADMIXTURE`
`ADMIXTURE` requires that you specify a `K` value (the number of bins). I ran `K = 1` to `K = 6` with 10 replicate runs per K and 10 cross-validation splits per run (a total of 60 runs with 10-fold CV for each).

```bash
# Run via job on UFRC, see admixture.job for details
# Resources used: 163 Mb, 10 min

module load admixture/1.23

export NAME="all_to_ASM1654582_fil_maf0.05"
# K = 1 through K = 6
for K in {1..6}
do
    # 10 replicate runs per K
    for r in {1..10}
    do
        echo doing K="$K",run="$r"
        admixture -s ${RANDOM} --cv=10 "$NAME".bed "$K" > log.K"$K".r"$r".out
        mv "$NAME"."$K".Q "$NAME".K"$K".r"$r".Q
        mv "$NAME"."$K".P "$NAME".K"$K".r"$r".P
    done
done
```
The results of each ADMIXTURE run can be viewed in the directory [admixture_output](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/admixture_output).

### `ADMIXTURE` post-processing

**Get cross-validation error values for each K:**
```bash
grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_to_ASM1654582_fil_maf0.05.cv.error
```

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
python vis_admixture.py
pong -m all_to_ASM1654582_fil_maf0.05_filemap.txt -i all_to_ASM1654582_fil_maf0.05_ind2pop.txt -n all_to_ASM1654582_fil_maf0.05_poporder.txt -l all_to_ASM1654582_fil_maf0.05_colors.txt
conda deactivate
```

See [`all_to_ASM1654582_fil_maf0.05_ADMIXTURE_pong`](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/all_to_ASM1654582_fil_maf0.05_ADMIXTURE_pong) for results of visualization.

## `ADMIXTURE` excluding sample WF03/1051
The individual WF03/1051 was distant from all other points in the PCA plots, so I ran ADMIXTURE without it to see whether that would change the results. I started by excluding this sample when generating the `BIM` file:

```bash
module load plink/1.90b3.39
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05"
PRUNED_VARS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca"

# Create file with samples to exclude
echo "WF03 WF03" > ids_to_exclude.txt

# Make BED and BIM files
plink --vcf "$VCF_DIR"/all_to_ASM1654582_fil_maf0.05.vcf --remove ids_to_exclude.txt --extract "$PRUNED_VARS_DIR"/all_maf0.05.prune.in --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out all_to_ASM1654582_fil_maf0.05
awk '{$1="0";print $0}' all_to_ASM1654582_fil_maf0.05.bim > all_to_ASM1654582_fil_maf0.05.bim.tmp
mv all_to_ASM1654582_fil_maf0.05.bim.tmp all_to_ASM1654582_fil_maf0.05.bim
```

All other steps were the same as listed above. You can find the files used in the subdirectory [without_wf03](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/without_wf03). There looked to be no difference between the run with and without WF03.