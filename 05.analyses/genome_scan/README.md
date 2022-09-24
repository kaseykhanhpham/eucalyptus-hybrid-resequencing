# Genome Scan Analysis

## Calculate nucleotide diversity statistics
First, I manually created lists of samples from the overall variant set that belonged to each population, under the names `Ecordata.txt`, `Eglobulus_MR.txt`, and `Eglobulus_ref.txt`. Each contains the RAPiD IDs under which each sample is identified in the VCF file, one per line, for each population. These will be used in all subsequent genome scan analyses to delimit populations.

Used [`vcftools`](https://vcftools.github.io) to calculate mean pairwise distance (pi) and Tajima's D in windows across the genome.

## Fst Calculations

**Calculate Fst in sliding windows along the genome:**

Calculated for Meehan Range _E. globulus_ and reference _E. globulus_ because I wanted to detect divergent regions in Meehan Range _E. globulus_ that may resemble _E. cordata_.

```bash
# Performed in UFRC queue system. See fst_windows.job for more details.
# Resources: 1.25 Mb, 1 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00_snps.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

vcftools --vcf "$INFILE" --out all_maf0.05 --weir-fst-pop "$POPLIST_DIR"/Eglobulus_MR.txt --weir-fst-pop "$POPLIST_DIR"/Eglobulus_ref.txt --fst-window-size 20000 --fst-window-step 2000
```

## Sliding Window Pi, Tajima's D, Dxy Calculations

Split multi-allelic SNPs into biallelic SNPs, remove problematic FORMAT fields which will cause errors:

```bash
# Performed in UFRC queue system. See format_egg.job for more details.
# Resources: 13 Mb, 2 min
module load bcftools/1.15

VCFDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00"

bcftools annotate -x FORMAT/PL,FORMAT/AD,FORMAT/AO,FORMAT/QA -O v -o - "$VCFDIR"/all_to_ASM1654582_fil_maf0.00_snps.vcf | bcftools norm -m -any -O v -o - > all_to_ASM1654582_fil_maf0.00_snps_biallelic_formatted.vcf
```

**Created list of chromosomes to iterate through**
Had to remove contigs which didn't have variants on them by hand before using this chromosome list file in the following script.

**Run python script to calculate Dxy in sliding windows across genome**
```bash
# Performed in UFRC queue system. See window_stats.job for more details.
# Resources: 210 Mb, 8 min
module load conda

ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

conda activate "$ENV_DIR"/euc_hyb_reseq

python "$SCRIPT_DIR"/stats_windows.py "$WDIR"/all_to_ASM1654582_fil_maf0.00_snps_biallelic_formatted.vcf 100000 20000 "$WDIR"/pop_structure.json "$WDIR"/outgroup_structure.json "$WDIR"/chr_list.txt
```

# Introgression Statistics

# DO QUIBL