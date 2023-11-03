# Coalescent Demographic History 

## MSMC (NOPE!!! NOPE!!)
Randomly picked 3 samples of Meehan Range _E. globulus_ and _E. cordata_ to use in MSMC analysis.
```R
# E. globulus
sample(1:20, 3, replace=FALSE)
# 11 1 4

# E. cordata
sample(1:10, 3, replace=FALSE)
# 3 2 10
```
Samples selected:
| Taxon            | RAPiD ID | Accession |
| ---------------- | -------- | --------- |
| _E. globulus_    | WA01     | 4190      |
| _E. globulus_    | WB02     | 6024      |
| _E. globulus_    | WE02     | 4224      |
| _E. cordata_     | WB01     | 5506      |
| _E. cordata_     | WD02     | 2899a     |
| _E. cordata_     | WH05     | 5509      |

Subsetted VCF file to these individuals.

## Meehan Range E. globulus
### Stairway Plot
Generated Site Frequency Spectrum using [easySFS](https://github.com/isaacovercast/easySFS). Manually generated population list file.

```bash
module load python/3.8
ESFS_DIR="/blue/soltis/kasey.pham/bin/easySFS"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -a -p esfs_poplist.txt --preview --total-length 8732510
# globulus segregating sites peaked (876054) at subsampling to 32 haplotypes
# cordata segregating sites peaked (690386) at subsampling to 16 haplotypes
"$ESFS_DIR"/easySFS.py -i "$VCF_DIR"/all_fil.vcf.gz -p esfs_poplist.txt -o esfs_outp -a -f --order glob_MR,cord_MR --proj 32,16 --total-length 8732510 -v
```

Created `Stairway Plot 2` input file and ran program.

```bash
module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder glob_fold.blueprint
bash glob_fold.blueprint.sh
```

## E. cordata
### Stairway Plot 2
Used Site Frequency Spectrum generated during _E. globulus_ analysis to make blueprint file.

```bash
module load java/20
PROG_DIR="/blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1"

java -cp "$PROG_DIR"/stairway_plot_es Stairbuilder cord_fold.blueprint
bash cord_fold.blueprint.sh > stairwayplot2.out
```
