# Genome Scan Analysis

## Calculate nucleotide diversity statistics
First, I manually created lists of samples from the overall variant set that belonged to each population, under the names `Ecordata.txt`, `Eglobulus_MR.txt`, and `Eglobulus_ref.txt`. Each contains the RAPiD IDs under which each sample is identified in the VCF file, one per line, for each population. These will be used in all subsequent genome scan analyses to delimit populations.

Used [`vcftools`](https://vcftools.github.io) to calculate mean pairwise distance (pi) and Tajima's D in windows across the genome.

```bash
# Performed in UFRC queue system. See pi.job for more details.
# Resources: 16 Mb, 2 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

vcftools --vcf "$INFILE" --out all_maf0.05 --window-pi 20000 --window-pi-step 2000 --keep "$POPLIST_DIR"/Eglobulus_MR.txt
```

```bash
# Performed in UFRC queue system. See tajimas_d.job for more details.
# Resources: 4.3 Mb, 2 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

vcftools --vcf "$INFILE" --out all_maf0.05 --TajimaD 10000 --keep "$POPLIST_DIR"/Eglobulus_MR.txt
```

## Divergence Calculations

### Fst Calculations

**Calculate Fst in sliding windows along the genome:**

Calculated for Meehan Range _E. globulus_ and reference _E. globulus_ because I wanted to detect divergent regions in Meehan Range _E. globulus_ that may resemble _E. cordata_.

```bash
# Performed in UFRC queue system. See fst_windows.job for more details.
# Resources: 4.3 Mb, 2 min

module load vcftools/0.1.16 

INFILE="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/filter_snps/maf0.00/all_to_ASM1654582_fil_maf0.00.vcf"
POPLIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"

vcftools --vcf "$INFILE" --out all_maf0.05 --weir-fst-pop "$POPLIST_DIR"/Eglobulus_MR.txt --weir-fst-pop "$POPLIST_DIR"/Eglobulus_ref.txt --fst-window-size 20000 --fst-window-step 2000
```

# Dxy Calculations

**Filter _E. cordata_ samples from MAF0.00 VCF**
```bash

```

Run python script to calculate Dxy in sliding windows across genome
```bash

```

# Introgression Statistics

# DO QUIBL