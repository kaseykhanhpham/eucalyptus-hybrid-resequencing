# Whole-genome `ADMIXTURE`
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/ADMIXTURE/).
Uses the maximum-likelihood population structure program [`ADMIXTURE`](http://dalexander.github.io/admixture/).

## `ADMIXTURE` on all ingroup samples

Done for both MAF=0.00 and MAF=0.05 variant sets.

### Create `BIM` file with variants for `ADMIXTURE`
This analysis needs a set of unlinked SNPs as input. I already generated a set of variants pruned for linkage for the PCA analysis, see [that section](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/PCA#prune-linked-snps) for more details.

```bash
module load plink/1.90b3.39
# maf=0.00
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.00"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf0.00.prune.in"
OUTNAME="meehan_all_fil_maf0.00"
INGROUP="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/ingroup.fam"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --keep "$INGROUP" --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim

# maf=0.05
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"
OUTNAME="meehan_all_fil_maf0.05"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --keep "$INGROUP" --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim
```

### Run `ADMIXTURE`
`ADMIXTURE` requires that you specify a `K` value (the number of bins). I ran `K = 1` to `K = 6` with 10 replicate runs per K and 10 cross-validation splits per run (a total of 60 runs with 10-fold CV for each). Showing the command for only MAF=0.05 below, for commands for each MAF, see `admixture_maf0.0X.job`.

```bash
# Run via job on UFRC, see admixture_maf0.05.job for details
# Resources used: 150 Mb, 1 hr (420 Mb, 3 hr for MAF=0.00)

module load admixture/1.23

export NAME="meehan_all_fil_maf0.05"
# K = 1 through K = 6
for K in {1..6}
do
    # 10 replicate runs per K
    for r in {1..10}
    do
        echo doing K="$K",run="$r"
        admixture -s ${RANDOM} --cv=10 "$NAME".bed "$K" > admixture_output/log.K"$K".r"$r".out
        mv "$NAME"."$K".Q admixture_output/"$NAME".K"$K".r"$r".Q
        mv "$NAME"."$K".P admixture_output/"$NAME".K"$K".r"$r".P
    done
done
```
The results of each ADMIXTURE run can be viewed in the directory admixture_output.
* [MAF=0.00](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/admixture_output)
* [MAF=0.05](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/admixture_output)

### `ADMIXTURE` post-processing

**Get cross-validation error values for each K:**
```bash
MAF00_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.00"
MAF05_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"

cd "$MAF00_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  meehan_all_fil_maf0.00.cv.error

cd "$MAF05_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  meehan_all_fil_maf0.05.cv.error
```

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00.py
pong -m meehan_all_fil_maf0.00_filemap.txt -i meehan_all_fil_maf0.00_ind2pop.txt -n meehan_all_fil_maf0.00_poporder.txt -l meehan_all_fil_maf0.00_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05.py
pong -m meehan_all_fil_maf0.05_filemap.txt -i meehan_all_fil_maf0.05_ind2pop.txt -n meehan_all_fil_maf0.05_poporder.txt -l meehan_all_fil_maf0.05_colors.txt

conda deactivate
```

Results of visualization:

![MAF=0.00 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), both groups of _E. globulus_ samples (blue) have estimated 0% admixture with _E. cordata_ (yellow). In K=3 (second lowest CV error), individuals across all _E. globulus_ sampling were assigned to one of two bins with little admixture inferred between. Suggests that there is probably a lot of variation to partition in _E. globulus_ as compared with _E. cordata_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/maf0.00_admixture.png "MAF = 0.00 ADMIXTURE K=2 through K=6")

![MAF=0.05 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), introgressed _E. globulus_ has an estimated 0% admixture rate with _E. cordata_ while "pure" _E. globulus_ has an estimated 2% admixture rate with _E. cordata_. I think this falls within negligible amounts attributable to being an ADMIXTURE artifact. In K=3 (second-lowest CV error), Meehan Range _E. globulus_ is inferred to be admixed between two bins, one of which all of reference _E. globulus_ is inferred to belong to.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/maf0.05_admixture.png "ADMIXTURE K=2 through K=6")

The K=3 result was surprising, as our reference panel for _E. globulus_ was trees from different populations and therefore probably higher diversity than those sampled from Meehan Range. However, because this pattern is not present in MAF=0.00, it can be inferred that it is being driven by the lack of singletons and increased number of shared variants. When filtering to shared variants, the more diverse reference pool would have rare or singleton alleles differentiating them from their brethren removed, leaving the dataset with the appearance that the reference pool is less diverse.

## `ADMIXTURE` excluding sample WF03/1051
The samples WF03/1051 was distant from all other points in the PCA plots, so I ran ADMIXTURE without it to see whether that would change the results. I started by excluding this sample when generating the `BIM` file:

```bash
module load plink/1.90b3.39
VCF00_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
VCF05_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"

PRUNED_VARS00_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf0.00.prune.in"
PRUNED_VARS05_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"

MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.05"

REMOVE_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/outlier_check/to_remove.fam"
# Used same file as in PCA analysis to indicate what samples to remove -- WF03/1051 and the outgroup.

# Make BED and BIM files
# MAF = 0.00
plink --vcf "$VCF00_FILE" --remove "$REMOVE_FILE" --extract "$PRUNED_VARS00_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$MAF00_DIR"/meehan_all_fil_maf0.00_outl
awk '{$1="0";print $0}' "$MAF00_DIR"/meehan_all_fil_maf0.00_outl.bim > "$MAF00_DIR"/meehan_all_fil_maf0.00_outl.bim.tmp
mv "$MAF00_DIR"/meehan_all_fil_maf0.00_outl.bim.tmp "$MAF00_DIR"/meehan_all_fil_maf0.00_outl.bim

# MAF = 0.05
plink --vcf "$VCF05_FILE" --remove "$REMOVE_FILE" --extract "$PRUNED_VARS05_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$MAF05_DIR"/meehan_all_fil_maf0.05_outl
awk '{$1="0";print $0}' "$MAF05_DIR"/meehan_all_fil_maf0.05_outl.bim > "$MAF05_DIR"/meehan_all_fil_maf0.05_outl.bim.tmp
mv "$MAF05_DIR"/meehan_all_fil_maf0.05_outl.bim.tmp "$MAF05_DIR"/meehan_all_fil_maf0.05_outl.bim
```

### Run `ADMIXTURE`

The following was done for both MAF=0.00 and MAF=0.05. Jobs below are examples. See respective job files for more details.

```bash
# Run via job on UFRC, see admixture_maf0.0X_outl.job for details
# Resources used: 900 Mb, 7 hours

module load admixture/1.23

export NAME="meehan_all_fil_maf0.05_outl"
# K = 1 through K = 6
for K in {1..6}
do
    # 10 replicate runs per K
    for r in {1..10}
    do
        echo doing K="$K",run="$r"
        admixture -s ${RANDOM} --cv=10 "$NAME".bed "$K" > admixture_output/log.K"$K".r"$r".out
        mv "$NAME"."$K".Q admixture_output/"$NAME".K"$K".r"$r".Q
        mv "$NAME"."$K".P admixture_output/"$NAME".K"$K".r"$r".P
    done
done
```

### `ADMIXTURE` post-processing

**Get cross-validation error values for each K:**
```bash
MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.05"

cd "$MAF00_DIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > meehan_all_fil_maf0.00_outl.cv.error

cd "$MAF05_DIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > meehan_all_fil_maf0.05_outl.cv.error
```

K = 2 had the lowest cross-validation error for both datasets.

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X_outl.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\outlier_check\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\outlier_check\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00_outl.py
pong -m meehan_all_fil_maf0.00_outl_filemap.txt -i meehan_all_fil_maf0.00_outl_ind2pop.txt -n meehan_all_fil_maf0.00_outl_poporder.txt -l meehan_all_fil_maf0.00_outl_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05_outl.py
pong -m meehan_all_fil_maf0.05_outl_filemap.txt -i meehan_all_fil_maf0.05_outl_ind2pop.txt -n meehan_all_fil_maf0.05_outl_poporder.txt -l meehan_all_fil_maf0.05_outl_colors.txt

conda deactivate
```
![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 for MAF=0.00, summarized over 10 individual runs in pong (without highly divergent reference sample.) Results resemble dataset with WF03 included. See K=3 plots for other run results.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.00/maf0.00_outl_admixture.png "ADMIXTURE K=2 through K=6 MAF=0.00 without WF03/1051")

![ADMIXTURE K=3 excluding samples WF03/1051 for MAF=0.00, summarized over 10 individual runs in pong (without highly divergent reference sample.) Results consistent with those where WF03 was included.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.00/maf0.00_outl_admixture_k3.png "ADMIXTURE K=3 MAF=0.00 without WF03/1051")

![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 for MAF=0.05, summarized over 10 individual runs in pong. Results are consistent with those where WF03 was included.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/weird_sample_check/maf0.05/maf0.05_outl_admixture.png "ADMIXTURE K=2 through K=6 MAF=0.05 without WF03/1051")
