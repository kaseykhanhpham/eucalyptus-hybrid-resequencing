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
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf0.00.prune.in"
OUTNAME="meehan_all_fil_maf0.00"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim

# maf=0.05
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"
OUTNAME="meehan_all_fil_maf0.05"

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
![MAF=0.05 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), introgressed _E. globulus_ has an estimated 1% admixture rate with _E. cordata_ while "pure" _E. globulus_ has an estimated 2% admixture rate with _E. cordata_. In K=3 (second-lowest CV error), Meehan Range _E. globulus_ has a portion of its genetic diversity attributed to a third bin.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/meehan_all_fil_maf0.05_ADMIXTURE.png "ADMIXTURE K=2 through K=6")

![MAF=0.00 ADMIXTURE K=3 for 10 individual runs in pong. Unlike in MAF = 0.05, variation is evenly attributed to two source bins across both reference and Meehan Range _E. globulus_.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/meehan_all_fil_maf0.00_ADMIXTURE.png "ADMIXTURE K=2 through K=6")

## `ADMIXTURE` excluding sample WF03/1051 and WG04/2025
The samples WF03/1051 and WG04/2025 were distant from all other points in the PCA plots, so I ran ADMIXTURE without them to see whether that would change the results. I started by excluding this sample when generating the `BIM` file:

```bash
module load plink/1.90b3.39
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
PRUNED_VARS_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf0.05.prune.in"
MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/weird_sample_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/weird_sample_check/maf0.05"

# Create file with samples to exclude
echo "WF03 WF03" > ids_to_exclude.txt
echo "WG04 WG04" >> ids_to_exclude.txt

# Make BED and BIM files
# MAF = 0.00
plink --vcf "$VCF_FILE" --remove ids_to_exclude.txt --extract "$PRUNED_VARS_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$MAF00_DIR"/meehan_all_fil_maf0.00_noWeird
awk '{$1="0";print $0}' "$MAF00_DIR"/meehan_all_fil_maf0.00_noWeird.bim > "$MAF00_DIR"/meehan_all_fil_maf0.00_noWeird.bim.tmp
mv "$MAF00_DIR"/meehan_all_fil_maf0.00_noWeird.bim.tmp "$MAF00_DIR"/meehan_all_fil_maf0.00_noWeird.bim

# MAF = 0.05
plink --vcf "$VCF_FILE" --remove ids_to_exclude.txt --extract "$PRUNED_VARS_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$MAF05_DIR"/meehan_all_fil_maf0.05_noWeird
awk '{$1="0";print $0}' "$MAF05_DIR"/meehan_all_fil_maf0.05_noWeird.bim > "$MAF05_DIR"/meehan_all_fil_maf0.05_noWeird.bim.tmp
mv "$MAF05_DIR"/meehan_all_fil_maf0.05_noWeird.bim.tmp "$MAF05_DIR"/meehan_all_fil_maf0.05_noWeird.bim
```

### Run `ADMIXTURE`

The following was done for both MAF=0.00 and MAF=0.05. Jobs below are examples. See respective job files for more details.

```bash
# Run via job on UFRC, see admixture_maf0.05_noWeird.job for details
# Resources used: 200 Mb, 2 hrs

module load admixture/1.23

export NAME="meehan_all_fil_maf0.05_noWeird"
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
MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/weird_sample_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/weird_sample_check/maf0.05"

cd "$MAF00_DIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > meehan_all_fil_maf0.00_noWeird.cv.error

cd "$MAF05_DIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > meehan_all_fil_maf0.05_noWeird.cv.error
```

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X_noWeird.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\weird_sample_check\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_ADMIXTURE\weird_sample_check\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00_noWeird.py
pong -m meehan_all_fil_maf0.00_noWeird_filemap.txt -i meehan_all_fil_maf0.00_noWeird_ind2pop.txt -n meehan_all_fil_maf0.00_noWeird_poporder.txt -l meehan_all_fil_maf0.00_noWeird_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05_noWeird.py
pong -m meehan_all_fil_maf0.05_noWeird_filemap.txt -i meehan_all_fil_maf0.05_noWeird_ind2pop.txt -n meehan_all_fil_maf0.05_noWeird_poporder.txt -l meehan_all_fil_maf0.05_noWeird_colors.txt

conda deactivate
```
![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 and WG04/2025 for MAF=0.00, summarized over 10 individual runs in pong. Without highly divergent samples. MAF=0.00 plot looks similar to MAF=0.05 plot, indicating that singletons in weird samples may have been driving high diversity and "homogenization" between introgressants and reference _E. globulus_?](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/weird_sample_check/maf0.00/meehan_all_fil_maf0.00_noWeird.png "ADMIXTURE K=2 through K=6 MAF=0.00 without weird samples")

![ADMIXTURE K=3 excluding samples WF03/1051 and WG04/2025 for MAF=0.00, summarized over 10 individual runs in pong. Some of these alternate results are interesting.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/weird_sample_check/maf0.00/meehan_all_fil_maf0.00_noWeird_K3.png "ADMIXTURE K=3 MAF=0.00 without weird samples")

![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 and WG04/2025 for MAF=0.05, summarized over 10 individual runs in pong. Results for K=2 and K=3 are fairly consistent with runs including WF03 and WG04](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/weird_sample_check/maf0.05/meehan_all_fil_maf0.05_noWeird.png "ADMIXTURE K=2 through K=6 MAF=0.05 without weird samples")
