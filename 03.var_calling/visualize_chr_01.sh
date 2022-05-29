#!/bin/bash

# Bash script for calculating statistics of variants
# Usage: visualize_chr.sh <VCF File> <name>
module load vcftools

vcftools --vcf "$1" --freq2 --max-alleles 2 --out "$2"
sed -i 's/{FREQ}/A1\tA2/g' "$2".frq
vcftools --vcf "$1" --depth --out "$2"
vcftools --vcf "$1" --site-mean-depth --out "$2"
vcftools --vcf "$1" --site-quality --out "$2"
vcftools --vcf "$1" --missing-indv --out "$2"
vcftools --vcf "$1" --missing-site --out "$2"
vcftools --vcf "$1" --het --out "$2"

