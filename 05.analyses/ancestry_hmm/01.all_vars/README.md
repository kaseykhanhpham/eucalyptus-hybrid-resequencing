# HMM Inference of Local Ancestry - All
Performed using the tool [`Ancestry_HMM`](https://github.com/russcd/Ancestry_HMM). Thanks to Shelley Sianta for suggesting this analysis.

## MAF = 0.00
### Prepare inputs
#### Genotype count input file
First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts.job for more details.
# Resources used: 5 Mb, 50 min

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out "$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### AncestryHMM input file
Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in.job for more details.
# Resources used: 4 Gb, 30 min

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/maf00/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt "$LIST_DIR"/coverage.txt sample_files.txt gen_dists.tab all_ahmm_in.tab
```

### Run the program
Ran `Ancestry_HMM` using generated genotype counts input files and results from [`ADMIXTURE`](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/wg_ADMIXTURE) to provide contributions from source populations.

Initial run parameters:
* -i: input of genotype counts among reference and introgressed samples
* -s: ploidy of each sample, all diploid
* -a: two source populations, 0 (_E. globulus_) contributed 99% of variation and 1 (E. cordata) contributed 1%
* -p: first ancestry pulse for glob background, num generations set above limit to indicate starting background, 99% of current admixed genomes 
* -p: second ancestry pulse for cordy background, start at 1700 generations (assuming pleistocene admixture and generation times of 10 years) and optimize, 1% of current admixed genomes but optimize estimate
* -g: genotype counts provided rather than read pileups
* -b: do 100 bootstraps of 10,000 SNP blocks
* -ne: effective population size of the admixed population -- I don't know this, so I'm not going to try to provide it.

```bash
# Run in UFRC queue system; see ancestryhmm.job for more details.
# Resources used: 11 Gb, 3 hrs

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i all_ahmm_in.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

### Analyze results

#### Summary plots
Plotted posteriors by chromosome for each admixed sample locally.
```bash
# Run in UFRC queue system; see plot_posteriors.job for more details.
# Resources used: 600 Mb, 10 min

module load R/4.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/plot_posteriors.r "$WDIR"/posteriors/"$NAME".posterior "goldenrod1,green3,deepskyblue4"
done < "$WDIR"/Eglobulus_MR.txt
```

Plotted heatmaps of presence/absence of inferred introgressed regions across chromosomes. See `maf00/heatmaps/p0.95` and `maf00/heatmaps/p0.90`.

#### Overlap with other outlier windows

Retrieved windows with a posterior probability of > 95% for homozygous or heterozygous _E. cordata ancestry_.

```bash
# Run in UFRC queue system; see get_cord_posterior.job for more details.
# Resources used: 500 Mb, 7 min

module load R/4.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/get_cord_posterior.r "$WDIR"/posteriors/"$NAME".posterior "$NAME"_cord 0.95
done < "$WDIR"/Eglobulus_MR.txt
```

Converted posterior tables to genome annotation BED files.

```bash
module load R/4.2

WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/post_to_bed.r "$NAME"_cord_het_0.95.tab bedfiles/"$NAME"_cord_het_0.95.bed_tmp
    Rscript "$SCRIPT_DIR"/post_to_bed.r "$NAME"_cord_hom_0.95.tab bedfiles/"$NAME"_cord_hom_0.95.bed_tmp
done < "$WDIR"/Eglobulus_MR.txt
```

Merged Ancestry_HMM loci together if <500 bp apart using `BEDTools`.

```bash
module load bedtools/2.30.0
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

while read NAME
do
    bedtools merge -d 500 -i "$NAME"_cord_het_0.95.bed_tmp > "$NAME"_cord_het_0.95.bed
    bedtools merge -d 500 -i "$NAME"_cord_hom_0.95.bed_tmp > "$NAME"_cord_hom_0.95.bed
done < "$LIST_DIR"/Eglobulus_MR.txt

rm *.bed_tmp
```

Examined overlap between `Ancestry_HMM`-inferred _cordata_ regions, regions of MR _globulus_ with low dxy, and regions with high fdM.

```bash
module load bedtools/2.30.0
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
DXY_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/pixy/individuals/outlier_files/bedfiles"
AHMM_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/maf00/cord_loci/bedfiles"
DSUITE_FILES="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan/dsuite"

while read NAME
do
    # AHMM with DXY
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p95.bed > "$NAME"_ahmm_het_dxy_p95.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p90.bed > "$NAME"_ahmm_het_dxy_p90.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p95.bed > "$NAME"_ahmm_hom_dxy_p95.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DXY_FILES"/cord_"$NAME"_dxy_outl_p90.bed > "$NAME"_ahmm_hom_dxy_p90.bed

    # AHMM with f-stats
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DSUITE_FILES"/fDm_40_20_outliers_p05.bed > "$NAME"_ahmm_het_fDm.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_het_0.95.bed -b "$DSUITE_FILES"/df_40_20_outliers_p05.bed > "$NAME"_ahmm_het_df.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DSUITE_FILES"/fDm_40_20_outliers_p05.bed > "$NAME"_ahmm_hom_fDm.bed
    bedtools intersect -a "$AHMM_FILES"/"$NAME"_cord_hom_0.95.bed -b "$DSUITE_FILES"/df_40_20_outliers_p05.bed > "$NAME"_ahmm_hom_df.bed

done < "$LIST_DIR"/Eglobulus_MR.txt
```

#### Statistics in Introgressed Windows
Delimited introgressed windows more precisely and checked SNP calling stats and genome scan stats for putative regions. See maf00/char_windows.

## MAF = 0.05
### Prepare inputs
#### Genotype count input file
First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts_maf05.job for more details.
# Resources used: 6 Mb, 2 hrs

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf - --out gt_counts/cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf -  --out gt_counts/glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf - --out gt_counts/"$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### AncestryHMM input file
Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in_maf05.job for more details.
# Resources used: 500 Mb, 7 min

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/maf05/gt_counts"

Rscript "$SCRIPT_DIR"/get_gen_dists.r "$COUNTS_DIR"/glob_ref.frq.count "$MAP_DIR"/1060_LH_F2_manual_copy.map gen_dists.tab
python "$SCRIPT_DIR"/make_ahmm_in.py ref_files.txt "$LIST_DIR"/coverage.txt sample_files.txt gen_dists.tab all_ahmm_in.tab
```

### Run the program
Ran `Ancestry_HMM` using generated genotype counts input files and results from [`ADMIXTURE`](https://github.com/kaseykhanhpham/eucalyptus-hybrid-resequencing/tree/main/05.analyses/wg_ADMIXTURE) to provide contributions from source populations.

Initial run parameters:
* -i: input of genotype counts among reference and introgressed samples
* -s: ploidy of each sample, all diploid
* -a: two source populations, 0 (_E. globulus_) contributed 99% of variation and 1 (E. cordata) contributed 1%
* -p: first ancestry pulse for glob background, num generations set above limit to indicate starting background, 99% of current admixed genomes 
* -p: second ancestry pulse for cordy background, start at 1700 generations (assuming pleistocene admixture and generation times of 10 years) and optimize, 1% of current admixed genomes but optimize estimate
* -g: genotype counts provided rather than read pileups
* -b: do 100 bootstraps of 10,000 SNP blocks
* -ne: effective population size of the admixed population -- I don't know this, so I'm not going to try to provide it.

```bash
# Run in UFRC queue system; see ancestryhmm_maf05.job for more details.
# Resources used: 2 Gb, 7 min

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i all_ahmm_in.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

### Analyze results

#### Summary Plots
Plotted posteriors by chromosome for each admixed sample locally.
```bash
# Run in UFRC queue system; see plot_posteriors_maf05.job for more details.
# Resources used:

module load R/4.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

while read NAME
do
    Rscript "$SCRIPT_DIR"/plot_posteriors.r "$WDIR"/maf05/posteriors/"$NAME".posterior "goldenrod1,green3,deepskyblue4"
done < "$WDIR"/Eglobulus_MR.txt
```

Plotted heatmaps of presence/absence of inferred introgressed regions across chromosomes. See maf05/heatmaps.

#### Statistics in Introgressed Windows
Delimited introgressed windows more precisely and checked SNP calling stats and genome scan stats for putative regions. See maf05/char_windows.