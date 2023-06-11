# Whole-Genome ADMIXTURE
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
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf00.prune.in"
OUTNAME="all_fil_maf0.00"
INGROUP="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/ingroup.fam"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --keep "$INGROUP" --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim

# maf=0.05
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf05.prune.in"
OUTNAME="all_fil_maf0.05"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --keep "$INGROUP" --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim
```

Ran `ADMIXTURE` for K=2 through K=6 with 10 replicate runs for cross-validation.

```bash
# Run via job on UFRC, see admixture_maf0.0X.job for details. Example from admixture_maf0.00.job below.
# Maximum resources used: 420 Mb, 4 hrs

module load admixture/1.23

export NAME="meehan_all_fil_maf0.00"
for K in {1..6}
do
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
MAF00_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.00"
MAF05_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf0.05"

cd "$MAF00_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  all_fil_maf0.00.cv.error

cd "$MAF05_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  all_fil_maf0.05.cv.error
```

K=2 has the lowest cross-validation error for both the MAF=0.00 and MAF=0.05 datasets.

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00.py
pong -m all_fil_maf0.00_filemap.txt -i all_fil_maf0.00_ind2pop.txt -n all_fil_maf0.00_poporder.txt -l all_fil_maf0.00_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05.py
pong -m all_fil_maf0.05_filemap.txt -i all_fil_maf0.05_ind2pop.txt -n all_fil_maf0.05_poporder.txt -l all_fil_maf0.05_colors.txt

conda deactivate
```

### Visualize `ADMIXTURE` Results

MAF = 0.00:

![MAF=0.00 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), both groups of _E. globulus_ samples (blue) have estimated 0% admixture with _E. cordata_ (yellow). In K=3 (second lowest CV error), individuals of reference _E. globulus_ were inferred to be non-admixed and to belong to one bin mainly. Meehan Range _E. globulus_ were inferred to be admixed between the reference bin and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/all_maf0.00_admixture.png "MAF = 0.00 ADMIXTURE K=2 through K=6")

Just K = 3:

![MAF=0.00 ADMIXTURE K=3 summarized over 10 individual runs in pong. Almost all runs display a similar pattern to the one displayed in the main ADMIXTURE visualization. Reference _E. globulus_ is inferred to be unadmixed, while Meehan Range _E. globulus_ is inferred to be admixed between the reference group and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/all_maf0.00_admixture_k3.png "MAF = 0.00 ADMIXTURE K=3")

MAF = 0.05:

![MAF=0.05 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), introgressed _E. globulus_ has an estimated 1% admixture rate with _E. cordata_ while "pure" _E. globulus_ has an estimated 2.6% admixture rate with _E. cordata_. In K=3 (second-lowest CV error), reference _E. globulus_ is inferred to be mostly unadmixed, while Meehan Range _E. globulus_ is inferred to be admixed between the reference group and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_maf0.05_admixture.png "MAF = 0.05 ADMIXTURE K=2 through K=6")

Just K = 3:

![MAF=0.05 ADMIXTURE K=3 summarized over 10 individual runs in pong. Almost all runs display a similar pattern to the one displayed in the main ADMIXTURE visualization.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_maf0.05_admixture_k3.png "MAF = 0.05 ADMIXTURE K=3")

**Discussion:**

The admixture inferred for reference _E. globulus_ with _E. cordata_ in K=2 for both datasets is quite small and probably falls within range of being an artifact of the analysis. Since the level of introgression expected is low, its signal may be largely swamped in a genome-wide `ADMIXTURE` plot.

The MAF=0.00 K=3 result was surprising, as our reference panel for _E. globulus_ was trees from different populations and therefore probably higher diversity than those sampled from Meehan Range. It's therefore strange that they would be inferred as less admixed than the Meehan Range samples. It's possible that my SNP filtering and thinning schema restricted to shared variants, excluding rare ones which would differentiate reference _E. globulus_ samples from each other.

(This pattern only showed up in MAF=0.05 when inferring structure with only variant sites identified in FreeBayes, which we assumed to be driven by excluding variation not shared across reference samples for MAF = 0.05.)

## `ADMIXTURE` excluding sample WF03/1051
The sample WF03/1051 was distant from all other points in the PCA plots, so I ran `ADMIXTURE` without it to see whether that would change the results. (Note: run on a previous dataset where linkage disequilibrium was pruned at r2=0.1 instead of r2=0.2. Results should be fairly consistent between the two, so I didn't re-do this analysis.) I started by excluding this sample when generating the `BIM` file:

```bash
module load plink/1.90b3.39
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"

PRUNED_VARS00_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00/all_maf00.prune.in"
PRUNED_VARS05_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.05/all_maf05.prune.in"

MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/outlier_check/maf0.05"

REMOVE_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/outlier_check/to_remove.fam"
# Used same file as in PCA analysis to indicate what samples to remove -- WF03/1051 and the outgroup.

# Make BED and BIM files
# MAF = 0.00
plink --vcf "$VCF_FILE" --remove "$REMOVE_FILE" --extract "$PRUNED_VARS00_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$MAF00_DIR"/all_fil_maf0.00_outl
awk '{$1="0";print $0}' "$MAF00_DIR"/all_fil_maf0.00_outl.bim > "$MAF00_DIR"/all_fil_maf0.00_outl.bim.tmp
mv "$MAF00_DIR"/all_fil_maf0.00_outl.bim.tmp "$MAF00_DIR"/all_fil_maf0.00_outl.bim

# MAF = 0.05
plink --vcf "$VCF_FILE" --remove "$REMOVE_FILE" --extract "$PRUNED_VARS05_FILE" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --maf 0.05 --out "$MAF05_DIR"/all_fil_maf0.05_outl
awk '{$1="0";print $0}' "$MAF05_DIR"/all_fil_maf0.05_outl.bim > "$MAF05_DIR"/all_fil_maf0.05_outl.bim.tmp
mv "$MAF05_DIR"/all_fil_maf0.05_outl.bim.tmp "$MAF05_DIR"/all_fil_maf0.05_outl.bim
```

### Run `ADMIXTURE`

Run `ADMIXTURE` for K=2 through K=6 with 10 replicate runs for cross-validation, excluding WF03/1051.

```bash
# Run via job on UFRC, see admixture_maf0.0X_outl.job for details. Example from admixture_maf0.00_outl.job below.
# Maximum resources used:

module load admixture/1.23

export NAME="all_fil_maf0.00_outl"
for K in {1..6}
do
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
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_fil_maf0.00_outl.cv.error

cd "$MAF05_DIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > all_fil_maf0.05_outl.cv.error
```

K = 2 had the lowest cross-validation error for both datasets.

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X_outl.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\outlier_check\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\outlier_check\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf0.00_outl.py
pong -m all_fil_maf0.00_outl_filemap.txt -i all_fil_maf0.00_outl_ind2pop.txt -n all_fil_maf0.00_outl_poporder.txt -l all_fil_maf0.00_outl_colors.txt

cd $MAF05_DIR
python vis_admixture_maf0.05_outl.py
pong -m all_fil_maf0.05_outl_filemap.txt -i all_fil_maf0.05_outl_ind2pop.txt -n all_fil_maf0.05_outl_poporder.txt -l all_fil_maf0.05_outl_colors.txt

conda deactivate
```

### Results of `ADMIXTURE` without WF03/1051

MAF = 0.00:

![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 for MAF=0.00, summarized over 10 individual runs in pong. Meehan Range _E. globulus_ is inferred to be admixed between two source populations while the reference remains more coherent without WF03/1051.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.00/all_fil_maf0.00_outl_admixture.png "ADMIXTURE K=2 through K=6 MAF=0.00 without WF03/1051")

Just K = 3:

![ADMIXTURE K=3 excluding samples WF03/1051 for MAF=0.00, summarized over 10 individual runs in pong (without highly divergent reference sample.) Reference _E. globulus_ inferred to be mostly from one source population as compared with analysis including WF03/1051.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.00/all_fil_maf0.00_outl_admixture_k3.png "ADMIXTURE K=3 MAF=0.00 without WF03/1051")

MAF = 0.05:

![ADMIXTURE K=2 through K=6 excluding samples WF03/1051 for MAF=0.05, summarized over 10 individual runs in pong. _Similar to the pattern seen in MAF=0.00, where reference _E. globulus_ is mostly assigned to one bin and Meehan Range _E. globulus_ is inferred to be admixed between the reference bin and another bin.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.05/all_fil_maf0.05_outl_admixture.png "ADMIXTURE K=2 through K=6 MAF=0.05 without WF03/1051")

Just K = 3:

![ADMIXTURE K=3 excluding samples WF03/1051 for MAF=0.05, summarized over 10 individual runs in pong (without highly divergent reference sample). Most runs infer admixed Meehan Range _E. globulus_, but one run infers _E. cordata_ samples to be grouped in two bins.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/outlier_check/maf0.05/all_fil_maf0.05_outl_admixture_k3.png "ADMIXTURE K=3 MAF=0.05 without WF03/1051")

Given that there are slight differences in results, but similar patterns recovered with and without WF03/1051, I will do the genome scan for pi/dXY/FST both with and without the strange sample.