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
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf00"
VCF_FILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf00/all_maf00.prune.in"
OUTNAME="all_fil_maf00"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim

# maf=0.05
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf05"
PRUNED_VARS="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf05/all_maf05.prune.in"
OUTNAME="all_fil_maf05"

cd "$WDIR"
plink --vcf "$VCF_FILE" --extract "$PRUNED_VARS" --set-missing-var-ids @:# --allow-extra-chr --vcf-half-call m --make-bed --out "$OUTNAME"
awk '{$1="0";print $0}' "$OUTNAME".bim > "$OUTNAME".bim.tmp
mv "$OUTNAME".bim.tmp "$OUTNAME".bim
```

Ran `ADMIXTURE` for K=2 through K=6 with 10 replicate runs for cross-validation.

```bash
# Run via job on UFRC, see admixture_maf0.0X.job for details. Example from admixture_maf0.00.job below.
# Maximum resources used: 1 Gb, 9 hrs

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
MAF00_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf00"
MAF05_WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_admixture/maf05"

cd "$MAF00_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  all_fil_maf00.cv.error

cd "$MAF05_WDIR" 
grep "CV" admixture_output/*.out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' >  all_fil_maf05.cv.error
```

K=2 has the lowest cross-validation error for both the MAF=0.00 and MAF=0.05 datasets.

**Visualize `ADMIXTURE` results:**

Done locally on personal computer. See `vis_admixture_maf0.0X.py` for processing steps to generate pong formatting files.

```bash
conda activate euc_hyb_reseq
$MAF00_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\maf0.00"
$MAF05_DIR = "C:\Users\Kasey\OneDrive - University of Florida\Grad School Documents\Projects\eucalyptus-hybrid-resequencing\05.analyses\wg_admixture\maf0.05"

cd $MAF00_DIR
python vis_admixture_maf00.py
pong -m all_fil_maf00_filemap.txt -i all_fil_maf00_ind2pop.txt -n all_fil_maf00_poporder.txt -l all_fil_maf00_colors.txt

cd $MAF05_DIR
python vis_admixture_maf05.py
pong -m all_fil_maf05_filemap.txt -i all_fil_maf05_ind2pop.txt -n all_fil_maf05_poporder.txt -l all_fil_maf05_colors.txt

conda deactivate
```

### Visualize `ADMIXTURE` Results

MAF = 0.00:

![MAF=0.00 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), both groups of _E. globulus_ samples (blue) have estimated 0% admixture with _E. cordata_ (yellow). In K=3 (second lowest CV error), individuals of all _E. globulus_ were inferred to be admixed between the reference bin and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/all_maf00_admixture.png "MAF = 0.00 ADMIXTURE K=2 through K=6")

Just K = 3:

![MAF=0.00 ADMIXTURE K=3 summarized over 10 individual runs in pong. Almost all runs display a similar pattern to the one displayed in the main ADMIXTURE visualization. All _E. globulus_ samples were inferred to be admixed between the reference group and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.00/all_maf00_admixture_k3.png "MAF = 0.00 ADMIXTURE K=3")

MAF = 0.05:

![MAF=0.05 ADMIXTURE K=2 through K=6 summarized over 10 individual runs in pong. In K=2 (lowest CV error), introgressed _E. globulus_ has an estimated 0.1% admixture rate with _E. cordata_ while "pure" _E. globulus_ has an estimated 1% admixture rate with _E. cordata_. In K=3 (second-lowest CV error), reference _E. globulus_ is inferred to be mostly unadmixed, while Meehan Range _E. globulus_ is inferred to be admixed between the reference group and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_maf05_admixture.png "MAF = 0.05 ADMIXTURE K=2 through K=6")

Just K = 3:

![MAF=0.05 ADMIXTURE K=3 summarized over 10 individual runs in pong. Almost all runs display a similar pattern to the one displayed in the main ADMIXTURE visualization, where reference _E. globulus_ was inferred to be non-admixed while Meehan Range _E. globulus_ was inferred to be admixed between the reference bin and another.](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/blob/main/05.analyses/wg_ADMIXTURE/maf0.05/all_maf05_admixture_k3.png "MAF = 0.05 ADMIXTURE K=3")

**Discussion:**

The admixture inferred for reference _E. globulus_ with _E. cordata_ in K=2 for both datasets is quite small and probably falls within range of being an artifact of the analysis. Since the level of introgression expected is low, its signal may be largely swamped in a genome-wide `ADMIXTURE` plot.

The MAF=0.05 K=3 result was surprising, as our reference panel for _E. globulus_ was trees from different populations and therefore probably higher diversity than those sampled from Meehan Range. It's therefore strange that they would be inferred as less admixed than the Meehan Range samples. We think that by excluding rare variants, the reference _E. globulus_ "loses" its diversity first before Meehan Range _E. globulus_ because the samples are not from a single population and much of the variation will be singletons.
