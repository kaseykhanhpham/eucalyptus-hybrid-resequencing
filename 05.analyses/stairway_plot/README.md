# [Stairway Plots](https://doi.org/10.1186/s13059-020-02196-9) to infer demographic change

## Prune VCF MAF = 0.00 by linkage disequilibrium
I already generated a set of unlinked SNPs for the [whole genome PCA](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/PCA) in `PLINK`. I used `PLINK` to generate a pruned VCF from this file.

```bash
# Run on UFRC's queue system, see vcf2sfs.job for more details on run parameters.
# Resources used: 5 Gb, 19 sec

module load plink/1.90b3.39 

INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00_snps.vcf"
PRUNEDDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_pca/maf0.00"

plink --vcf "$INFILE00" --double-id --allow-extra-chr --set-missing-var-ids @:# --extract "$PRUNEDDIR"/all_maf0.00.prune.in --vcf-half-call m --recode vcf --out all_maf0.00
```

## Convert VCF file into Site Frequency Spectra
Used [Liu et al. 2018](https://doi.org/10.1111/mec.14782)'s [code](https://github.com/zhongmicai/Vcf2SiteFrequencySpectrum) for converting VCF to SFS.

For the population map file `glob-cord_popmap.txt`, used the following integer IDs for each population:

| population      | id |
| --------------- | -- |
| E. globulus MR  | 1  |
| E. globulus ref | 2  |
| E. cordata ref  | 3  |

```R
# Run on UFRC's queue system, see vcf2sfs.job for more details on run parameters.
# Resources used: 1.59 Mb, 30 sec

script_dir <- "/blue/soltis/kasey.pham/bin"
w_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/stairway_plot"
source(paste(script_dir, "vcf2sfs.r", sep = "/"))

setwd(w_dir)

# Run for each population
# E. globulus introgressants
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "glob_MR.fs", 1)

# E. globulus reference
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "glob_pure.fs", 2)

# E. cordata reference
vcf2dadi("all_maf0.00_pruned.vcf", "glob-cord_popmap.txt", "cord.fs", 3)
```

## Generate Stairway Plots
Used [Liu and Fu 2020's](https://doi.org/10.1186%2Fs13059-020-02196-9) [`Stairway Plot 2` program](https://github.com/xiaoming-liu/stairway-plot-v2).

Mutation rate for _E. grandis_ estimated in [Silva-Junior and Grattapaglia 2015](https://doi.org/10.1111/nph.13505) used in blueprint files.

Assuming a 10 year generation time based on the one given for wild _E. grandis_ in Eldridge et al. 1993, but will check a physical copy at the library and revise if incorrect.

Get information for blueprint files:

```bash
# number of variants
module load bcftools 
bcftools stats  all_maf0.00_pruned.vcf # 276063
```

Generate batch file to run the program:
```bash
module load java/14

java -cp /blue/soltis/kasey.pham/bin/stairway_plot_v2.1.1/stairway_plot_es Stairbuilder glob_mr.blueprint
```