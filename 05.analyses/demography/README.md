# Coalescent Demographic History 

## Stairway Plot
### _E. globulus_
Generated Site Frequency Spectrum using [easySFS](https://github.com/isaacovercast/easySFS). Manually generated population list file.

```bash
module load python/3.8
ESFS_DIR="/blue/soltis/kasey.pham/bin/easySFS"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -a -p esfs_poplist.txt --preview
# globulus segregating sites peaked (876054) at subsampling to 32 haplotypes
# cordata segregating sites peaked (690386) at subsampling to 16 haplotypes
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -p esfs_poplist.txt -o esfs_outp -a -f --order glob_MR,cord_MR --proj 32,16 -v
```

Created `Stairway Plot 2` input file and ran program.

```bash
# Ran in UFRC queue system; see stairwayplot2_glob.job for more details.
# Resources used:

module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder glob_fold.blueprint
bash glob_fold.blueprint.sh
```

### _E. cordata_
Used Site Frequency Spectrum generated during _E. globulus_ analysis to make blueprint file.

```bash
# Ran in UFRC queue system; see stairwayplot2_cord.job for more details.
# Resources used:
module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder cord_fold.blueprint
bash cord_fold.blueprint.sh
```

## Demography model fitting
Used [`dadi`](https://dadi.readthedocs.io/en/latest/) to fit different demography models to observed variant distribution.
Models fitted:
* Divergence, no migration, with size change
* Divergence with asymmetric migration, isolation, with size change
* Divergence in isolation, secondary contact, with size change
* Divergence in isolation, secondary contact, isolation, with size change

First created folded SFS for observed data and bootstrapped variants in `dadi`. Followed the manual's [2D demography example](https://dadi.readthedocs.io/en/latest/examples/basic_workflow/basic_workflow_2d_demographics/).
```bash
module load dadi
python dadi_prep.py
```

Ran parameter optimization for all 2D demographic models as described in Portik et al. 2017 and the associated pipeline, [`dadi_pipeline`](https://github.com/dportik/dadi_pipeline/tree/master). See individual python files and `dadi_pipeline`'s documentation for more details. In short, ran suggested 4 round optimization, using the previous round's best replicate (judged by likelihood) as a seed for heuristic search starting points with decreasing amounts of perturbation.