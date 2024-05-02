# HMM Inference of Local Ancestry
Performed using the tool [`Ancestry_HMM`](https://github.com/russcd/Ancestry_HMM). Thanks to Shelley Sianta for suggesting this analysis.

## MAF = 0.00
### Prepare inputs
#### Genotype count input file
First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts_maf00.job for more details.
# Resources used: 10 Mb, 1.5 hrs

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --out gt_counts/"$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### AncestryHMM input file
Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in_maf00.job for more details.
# Resources used: 5.8 Gb, 1.5 hrs

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/01.all_vars/maf00/gt_counts"

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
# Run in UFRC queue system; see ancestryhmm_maf00.job for more details.
# Resources used: 14 Gb, 2 hrs

module load ancestryhmm/1.0.2
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i all_ahmm_in.tab -s "$LIST_DIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

### Analyze results
#### Heatmap
See `heatmap` directory README for code on creating summary heatmaps.

#### Stats of ROIs
Calculated average pi, dxy, and FST in AHMM ROIs.

## MAF = 0.05
### Prepare inputs
#### Genotype count input file
First retrieved allele counts for each sample using `VCFTools`. Excluded ChrUn from this analysis as it does not have corresponding distances in the genetic map.

```bash
# Run in UFRC queue system; see gt_counts_maf05.job for more details.
# Resources used: 13 Mb, 2.5 hrs

module load vcftools/0.1.16 

IN_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

# Get genotype counts for E. cordata
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf - --out gt_counts/cord_ref --keep "$POPLIST_DIR"/Ecordata.txt --min-alleles 1 --max-alleles 2 --maf 0.05 --not-chr ChrUn --counts

# E. globulus ref
vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf -  --out gt_counts/glob_ref --keep "$POPLIST_DIR"/Eglobulus_ref.txt --min-alleles 1 --max-alleles 2 --maf 0.05 --not-chr ChrUn --counts

# Introgressants
while read NAME
do
    vcftools --gzvcf "$IN_DIR"/all_fil.vcf.gz --maf 0.05 --recode --stdout | vcftools --vcf - --out gt_counts/"$NAME" --indv "$NAME" --min-alleles 1 --max-alleles 2 --maf 0.05 --not-chr ChrUn --counts
done < "$POPLIST_DIR"/Eglobulus_MR.txt
```

### AncestryHMM input file
Generated input file for `Ancestry_HMM` using a custom `python` script.

```bash
# Run in UFRC queue system; see get_ahmm_in05.job for more details.
# Resources used: 630 Mb, 8 min

module load R/4.2

LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/flare"
COUNTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm/01.all_vars/maf05/gt_counts"

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
# Resources used: 1.6 Gb, 10 min

module load ancestryhmm/1.0.2
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ancestry_hmm"

ancestry_hmm -i "$WDIR"/all_ahmm_in.tab -s "$WDIR"/sample_ploidy.txt -a 2 0.99 0.01 -p 0 100000 0.99 -p 1 -1700 0.01 -g -b 100 1000
```

## Characterized windows by statistics

See maf00/char_windows and maf05/char_windows.