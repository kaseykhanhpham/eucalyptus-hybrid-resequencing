# Linkage Disequilibrium
## Allele Phasing
Index VCF files
```bash
module load samtools/0.1.18
ls *.bam | while read FILE
do
    tabix "$FILE"
done
```

Ran WhatsHap
```bash
# Ran in UFRC queue system; see whatshap.job for more details.
# Resources Used:

module load whatshap/1.1
REF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
BAM_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/02.process_reads/04.markdup"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps"
declare -a VCFLIST=(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chrUn)

for NAME in "${VCFLIST[@]}"
do
    whatshap phase -o "$NAME"_hapmarked.vcf --reference "$REF_DIR"/EGLOB-X46.v1.0.fa --tag PS "$VCF_DIR"/"$NAME"_fil.vcf.gz "$BAM_DIR"/S10_marked.bam "$BAM_DIR"/S2_marked.bam "$BAM_DIR"/S337_marked.bam "$BAM_DIR"/S356_marked.bam "$BAM_DIR"/S11_marked.bam "$BAM_DIR"/S30_marked.bam "$BAM_DIR"/S338_marked.bam "$BAM_DIR"/S357_marked.bam "$BAM_DIR"/S12_marked.bam "$BAM_DIR"/S31_marked.bam "$BAM_DIR"/S339_marked.bam "$BAM_DIR"/S358_marked.bam "$BAM_DIR"/S13_marked.bam "$BAM_DIR"/S320_marked.bam "$BAM_DIR"/S33_marked.bam "$BAM_DIR"/S359_marked.bam "$BAM_DIR"/S14_marked.bam "$BAM_DIR"/S321_marked.bam "$BAM_DIR"/S340_marked.bam "$BAM_DIR"/S35_marked.bam "$BAM_DIR"/S15_marked.bam "$BAM_DIR"/S322_marked.bam "$BAM_DIR"/S341_marked.bam "$BAM_DIR"/S36_marked.bam "$BAM_DIR"/S16_marked.bam "$BAM_DIR"/S323_marked.bam "$BAM_DIR"/S342_marked.bam "$BAM_DIR"/S37_marked.bam "$BAM_DIR"/S17_marked.bam "$BAM_DIR"/S324_marked.bam "$BAM_DIR"/S343_marked.bam "$BAM_DIR"/S38_marked.bam "$BAM_DIR"/S18_marked.bam "$BAM_DIR"/S325_marked.bam "$BAM_DIR"/S344_marked.bam "$BAM_DIR"/S39_marked.bam "$BAM_DIR"/S19_marked.bam "$BAM_DIR"/S326_marked.bam "$BAM_DIR"/S345_marked.bam "$BAM_DIR"/S3_marked.bam "$BAM_DIR"/S1_marked.bam "$BAM_DIR"/S327_marked.bam "$BAM_DIR"/S346_marked.bam "$BAM_DIR"/S40_marked.bam "$BAM_DIR"/S20_marked.bam "$BAM_DIR"/S328_marked.bam "$BAM_DIR"/S347_marked.bam "$BAM_DIR"/S4_marked.bam "$BAM_DIR"/S21_marked.bam "$BAM_DIR"/S329_marked.bam "$BAM_DIR"/S348_marked.bam "$BAM_DIR"/S5_marked.bam "$BAM_DIR"/S22_marked.bam "$BAM_DIR"/S32_marked.bam "$BAM_DIR"/S349_marked.bam "$BAM_DIR"/S6_marked.bam "$BAM_DIR"/S23_marked.bam "$BAM_DIR"/S330_marked.bam "$BAM_DIR"/S34_marked.bam "$BAM_DIR"/S7_marked.bam "$BAM_DIR"/S24_marked.bam "$BAM_DIR"/S331_marked.bam "$BAM_DIR"/S350_marked.bam "$BAM_DIR"/S8_marked.bam "$BAM_DIR"/S25_marked.bam "$BAM_DIR"/S332_marked.bam "$BAM_DIR"/S351_marked.bam "$BAM_DIR"/S9_marked.bam "$BAM_DIR"/S26_marked.bam "$BAM_DIR"/S333_marked.bam "$BAM_DIR"/S352_marked.bam "$BAM_DIR"/SRR10339635_marked.bam "$BAM_DIR"/S27_marked.bam "$BAM_DIR"/S334_marked.bam "$BAM_DIR"/S353_marked.bam "$BAM_DIR"/S28_marked.bam "$BAM_DIR"/S335_marked.bam "$BAM_DIR"/S354_marked.bam "$BAM_DIR"/S29_marked.bam "$BAM_DIR"/S336_marked.bam "$BAM_DIR"/S355_marked.bam
done
```

Compressed and indexed WhatsHap output.
```bash
module load bcftools
for NAME in "${VCFLIST[@]}"; do bgzip "$NAME"_hapmarked.vcf; bcftools index "$NAME"_hapmarked.vcf.gz; done
```

Reformatted genetic map into HapMap format and then ran SHAPEIT4.
```bash
# Ran in UFRC queue system; see shapeit.job for more details.
# Resources Used: 330 Mb, 4 min

module load shapeit4/4.2.2
MAP_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
declare -a VCFLIST=(01 02 03 04 05 06 07 08 09 10 11)

for NUM in "${VCFLIST[@]}"
do
    shapeit4 --input chr"$NUM"_hapmarked.vcf.gz --map "$MAP_DIR"/1060_LH_F2_manual_copy_Chr"$NUM".gmap --region Chr"$NUM" --output chr"$NUM"_phased.vcf --thread 12
done
```
## E. globulus LD
Used [`emeraLD`](https://github.com/statgen/emeraLD/tree/master) to calculate LD in 100kb non-overlapping sliding windows across all chromosomes.

Used `tabix` to index phased variants
```bash
module load htslib/1.15
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
declare -a VCFLIST=(01 02 03 04 05 06 07 08 09 10 11)

cd "$WDIR"
for NUM in "${VCFLIST[@]}"
do
    bgzip chr"$NUM"_hapmarked.vcf
    bgzip chr"$NUM"_phased.vcf
    tabix chr"$NUM"_phased.vcf.gz
done
```

Ran `emeraLD` across each chromosome for all pairwise SNPs less than 1Mbp apart.
```bash
# Ran on UFRC queue system; see genwide_ld_mac01.job for more details.
# Resources used:


```

Generated jobs to run emeraLD across windows in each chromosome
```bash
module load python/3.8
SCRIPT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
VCF_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/phase"
LIST_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/ld/glob"

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr01_phased.vcf.gz Chr01 42219553 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr01/mac01/chr01_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr02_phased.vcf.gz Chr02 50828380 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr02/mac01/chr02_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr03_phased.vcf.gz Chr03 65547241 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr03/mac01/chr03_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr04_phased.vcf.gz Chr04 38599333 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr04/mac01/chr04_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr05_phased.vcf.gz Chr05 62919971 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr05/mac01/chr05_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr06_phased.vcf.gz Chr06 52140660 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr06/mac01/chr06_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr07_phased.vcf.gz Chr07 54252628 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr07/mac01/chr07_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr08_phased.vcf.gz Chr08 70214608 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr08/mac01/chr08_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr09_phased.vcf.gz Chr09 38300324 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr09/mac01/chr09_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr10_phased.vcf.gz Chr10 38722660 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr10/mac01/chr10_ld_mac01.job

python "$SCRIPT_DIR"/make_emerald_job.py "$VCF_DIR"/chr11_phased.vcf.gz Chr11 42056460 100000 1 "$LIST_DIR"/Eglobulus_MR.txt "$LIST_DIR"/emerald/chr11/mac01/chr11_ld_mac01.job
```

Not displaying example of jobfiles here since they are quite large; see chrXX_ld_mac01.job in the jobfiles directory.

## E. cordata LD
