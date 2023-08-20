# Local Ancestry Inference with `ELAI`
Used the local ancestry inference software [`Efficient Local Ancestry Inference`](https://github.com/haplotype/elai) to infer admixture breakpoints along genome. This software uses a two-level HMM model based on the multispecies coalescent, with the upper layer corresponding to different source populations and the lower level corresponding to haplotype variation within a population.

## Input File Preparation
Used `PLINK` to convert VCF files to `BIMBAM` format genotype files using the [`recode`](https://www.cog-genomics.org/plink/1.9/data#recode) function.

```bash
# Run on UFRC queue system; see make_bimbam.job for more details.
# Resources used: 2 Gb, 2 min

module load plink/1.90b3.39 

VCF="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/all_fil.vcf.gz"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"

# MR cordata
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Ecordata.fam --out "$WDIR"/cord --recode bimbam

# ref globulus
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Eglobulus_ref.fam --out "$WDIR"/glob_pure --recode bimbam

# MR globulus
plink --vcf "$VCF" --double-id --allow-extra-chr 0 --set-missing-var-ids @:# --keep "$WDIR"/Eglobulus_MR.fam --out "$WDIR"/glob_mr --recode bimbam
```

Split the genotype files by chromosome (as per the recommendation of the `ELAI` manual).

```bash
# Run on UFRC queue system; see split_bimbam.job for more details.
# Resources used: 17 Mb, 2 sec

module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
python "$SCRIPTS_DIR"/split_bimbam.py
```

Replaced the number of SNPs in headings

```bash
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai/geno_by_chr"
cd "$WDIR"

declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

wc -l pos_Chr*.txt
declare -a SIZELIST=(787645 839980 736729 620251 620838 1079716 661655 1037845 688961 780227 774355)

rename _geno.txt _geno.tmp *

for INDEX in {0..10};
do
    head -n 1 cord_"${CHRLIST[$INDEX]}"_geno.tmp > cord_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> cord_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 cord_"${CHRLIST[$INDEX]}"_geno.tmp >> cord_"${CHRLIST[$INDEX]}"_geno.txt

    head -n 1 glob_ref_"${CHRLIST[$INDEX]}"_geno.tmp > glob_ref_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> glob_ref_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 glob_ref_"${CHRLIST[$INDEX]}"_geno.tmp >> glob_ref_"${CHRLIST[$INDEX]}"_geno.txt

    head -n 1 glob_mr_"${CHRLIST[$INDEX]}"_geno.tmp > glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
    echo "${SIZELIST[$INDEX]}" >> glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
    tail -n +3 glob_mr_"${CHRLIST[$INDEX]}"_geno.tmp >> glob_mr_"${CHRLIST[$INDEX]}"_geno.txt
done

rm *.tmp
```

## Ran `ELAI`

For MAF = 0.00 and estimated generation time since admixture = 100 (based on discussion with Brad and Rene, re: 40k year chloroplast capture). Multiple runs performed and results averaged as per recommendations of the developer; example given below.

```bash
# Run on UFRC queue system; see elai_maf00_g100_r1.job for more details.
# Resources used:

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHRLIST[@]}"
do
    "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r1 -C 2 -c 10 -mg 100
done
```

For MAF = 0.05 and estimated generation time since admixture = 100. Multiple runs performed and results averaged as per recommendations of the developer; example given below.

```bash
# Run on UFRC queue system; see elai_maf05_g100_r1.job for more details.
# Resources used:

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHRLIST[@]}"
do
    "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r1 -C 2 -c 10 -mg 100 -exclude-maf 0.05
done
```

For MAF = 0.00 and estimated generation time since admixture = 800 (based on initial estimate given by `AncestryHMM`). Multiple runs performed and results averaged as per recommendations of the developer; example given below.

```bash
# Run on UFRC queue system; see elai_maf00_g800_r1.job for more details.
# Resources used:

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHRLIST[@]}"
do
    "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r1 -C 2 -c 10 -mg 800
done
```

For MAF = 0.05 and estimated generation time since admixture = 800. Multiple runs performed and results averaged as per recommendations of the developer; example given below.

```bash
# Run on UFRC queue system; see elai_maf05_g800_r1.job for more details.
# Resources used:

module load gsl/2.6

BIN_DIR="/blue/soltis/kasey.pham/bin/ELAI"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/elai"
declare -a CHRLIST=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

for CHR in "${CHRLIST[@]}"
do
    "$BIN_DIR"/elai -g "$WDIR"/geno_by_chr/glob_ref_"$CHR"_geno.txt -p 10 -g "$WDIR"/geno_by_chr/cord_"$CHR"_geno.txt -p 11 -g "$WDIR"/geno_by_chr/glob_mr_"$CHR"_geno.txt -p 1 -pos "$WDIR"/geno_by_chr/pos_"$CHR".txt -s 30 -o "$CHR"_r1 -C 2 -c 10 -mg 800 -exclude-maf 0.05
done
```