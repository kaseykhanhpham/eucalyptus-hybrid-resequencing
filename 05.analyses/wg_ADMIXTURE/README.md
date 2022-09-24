# Whole-genome `ADMIXTURE`
Procedure based heavily on protocol in [Ravinet and Meier's tutorial](https://speciationgenomics.github.io/ADMIXTURE/).
Uses the maximum-likelihood population structure program [`ADMIXTURE`](http://dalexander.github.io/admixture/).

## `ADMIXTURE` on all samples

Done for both MAF=0.00 and MAF=0.05 variant sets.

### Create `BIM` file with variants for `ADMIXTURE`
This analysis needs a set of unlinked SNPs as input. I already generated a set of variants pruned for linkage for the PCA analysis, see [that section](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/PCA#prune-linked-snps) for more details.

```bash
module load plink/1.90b3.39
# maf=0.00
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.00"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf0.00.prune.in"
OUTNAME="all_to_ASM1654582_fil_maf0.00"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim

# maf=0.05
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"
OUTNAME="all_to_ASM1654582_fil_maf0.05"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim
```

### Run `ADMIXTURE`
`ADMIXTURE` requires that you specify a `K` value (the number of bins). I ran `K = 1` to `K = 6` with 10 replicate runs per K and 10 cross-validation splits per run (a total of 60 runs with 10-fold CV for each). Showing the command for only MAF=0.05 below, for commands for each MAF, see `admixture_maf0.0X.job`.

```bash
# Run via job on UFRC, see admixture_maf0.05.job for details
# Resources used: 150 Mb, 1 hr (420 Mb, 3 hr for MAF=0.00)

module load admixture/1.23

export NAME="all_to_ASM1654582_fil_maf0.05"
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
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_to_ASM1654582_fil_maf0.00.cv.error

cd "$MAF05_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_to_ASM1654582_fil_maf0.05.cv.error
```

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00.py
pong -m all_to_ASM1654582_fil_maf0.00_filemap.txt -i all_to_ASM1654582_fil_maf0.00_ind2pop.txt -n all_to_ASM1654582_fil_maf0.00_poporder.txt -l all_to_ASM1654582_fil_maf0.00_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05.py
pong -m all_to_ASM1654582_fil_maf0.05_filemap.txt -i all_to_ASM1654582_fil_maf0.05_ind2pop.txt -n all_to_ASM1654582_fil_maf0.05_poporder.txt -l all_to_ASM1654582_fil_maf0.05_colors.txt

conda deactivate
```

Results of visualization:
![MAF=0.05 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), introgressed _E. globulus_ has an estimated 1% admixture rate with _E. cordata_ while "pure" _E. globulus_ has an estimated 2% admixture rate with _E. cordata_. In K=3 (second-lowest CV error), Meehan Range _E. globulus_ has a portion of its genetic diversity attributed to a third bin.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_to_ASM1654582_fil_maf0.05_pong.png "ADMIXTURE K=2 through K=6")

![MAF=0.00 ADMIXTURE K=3 for 10 individual runs in pong. Unlike in MAF = 0.05, variation is evenly attributed to two source bins across both reference and Meehan Range _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_to_ASM1654582_fil_maf0.05_pong.png "ADMIXTURE K=2 through K=6")

## `ADMIXTURE` excluding sample WF03/1051
The individual WF03/1051 was distant from all other points in the PCA plots, so I ran ADMIXTURE without it to see whether that would change the results. I started by excluding this sample when generating the `BIM` file:

```bash
module load plink/1.90b3.39
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.05/all_to_ASM1654582_fil_maf0.05_snps.vcf"
PRUNED_VARS_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"

# Create file with samples to exclude
echo "WF03 WF03" > ids_to_exclude.txt

# Make BED and BIM files
plink --vcf "$VCF_FILE" --remove ids_to_exclude.txt --extract "$PRUNED_VARS_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out all_to_ASM1654582_fil_maf0.05_noWF03
awk '{$1="0";print $0}' all_to_ASM1654582_fil_maf0.05_noWF03.bim > all_to_ASM1654582_fil_maf0.05_noWF03.bim.tmp
mv all_to_ASM1654582_fil_maf0.05_noWF03.bim.tmp all_to_ASM1654582_fil_maf0.05_noWF03.bim
```

### Run `ADMIXTURE`

```bash
# Run via job on UFRC, see admixture_maf0.05_noWF03.job for details
# Resources used: 150 Mb, 1 hr

module load admixture/1.23

export NAME="all_to_ASM1654582_fil_maf0.05_noWF03"
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
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05/without_WF03"

cd "$WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_to_ASM1654582_fil_maf0.05_noWF03.cv.error
```

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.05_noWF03.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$WDIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\without_wf03"

cd $WDIR
python vis_admixture_maf0.05_noWF03.py
pong -m all_to_ASM1654582_fil_maf0.05_noWF03_filemap.txt -i all_to_ASM1654582_fil_maf0.05_noWF03_ind2pop.txt -n all_to_ASM1654582_fil_maf0.05_noWF03_poporder.txt -l all_to_ASM1654582_fil_maf0.05_noWF03_colors.txt

conda deactivate
```

![ADMIXTURE K=2 through K=6 excluding sample WF03, summarized over 10 individual runs in pong. In K=2, introgressed _E. globulus_ still has an estimated 1% admixture rate with _E. cordata_ while "pure" _E. globulus_ still has an estimated 2% admixture rate with _E. cordata_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/without_wf03/all_to_ASM1654582_fil_maf0.05_noWF03_pong.png "ADMIXTURE K=2 through K=6 without WF03")
