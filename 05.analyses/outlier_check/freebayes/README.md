# Outlier sample check

## Wier-Cockerham FST
Checked values of FST between populations with and without strange samples which stood out in whole genome PCA to try to get a better understanding of the distribution of variation in samples. WF03/1051 was distant in 2nd and 3rd PC plots from the rest of the _E. globulus_ samples, so did separate calculations against the different populations and ommitting it in different combinations.

```bash
module load vcftools/0.1.16

INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
LISTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check"
MAF00_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check/maf0.00"
MAF05_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check/maf0.05"

# ALL SAMPLE GROUPS AGAINST EACH OTHER
# All, MAF=0.00
cd "$MAF00_DIR"
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_all_to_all
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.00_cord_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_cord_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_globmr_globref

# All, MAF=0.05
cd "$MAF05_DIR"
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_all_to_all
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.05_cord_globmr
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_cord_globref
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_globmr_globref

# SAMPLE GROUPS AGAINST EACH OTHER WITH PROBLEMATIC SAMPLES REMOVED
# MAF = 0.00
cd "$MAF00_DIR"
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_all_to_all_outl
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.00_cord_globmr_outl
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_cord_globref_outl
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_globmr_globref_outl

# MAF = 0.05
cd "$MAF05_DIR"
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_all_to_all_outl
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.05_cord_globmr_outl
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_cord_globref_outl
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_globmr_globref_outl
```

Get genome-wide average in `R`:

```R
maf00_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check/maf0.00"
maf05_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check/maf0.05"

# ALL SAMPLE GROUPS 
# MAF=0.00, all sample groups against each other
setwd(maf00_dir)
# import FST results
maf00_all_to_all <- read.table("maf0.00_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
maf00_cord_globmr <- read.table("maf0.00_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
maf00_cord_globref <- read.table("maf0.00_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
maf00_globmr_globref <- read.table("maf0.00_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")
# get summaries of FST distrs
summary(na.omit(maf00_all_to_all$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_cord_globmr$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_cord_globref$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_globmr_globref$WEIR_AND_COCKERHAM_FST))

# MAF=0.00, sample groups with WF03 removed against each other
# import FST results
maf00_all_to_all_outl <- read.table("maf0.00_all_to_all_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf00_cord_globmr_outl <- read.table("maf0.00_cord_globmr_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf00_cord_globref_outl <- read.table("maf0.00_cord_globref_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf00_globmr_globref_outl <- read.table("maf0.00_globmr_globref_outl.weir.fst", header = TRUE, na.strings = "-nan")
# get summaries of FST distrs
summary(na.omit(maf00_all_to_all_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_cord_globmr_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_cord_globref_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf00_globmr_globref_outl$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, all sample groups against each other
setwd(maf05_dir)
# import FST results
maf05_all_to_all <- read.table("maf0.05_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
maf05_cord_globmr <- read.table("maf0.05_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
maf05_cord_globref <- read.table("maf0.05_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
maf05_globmr_globref <- read.table("maf0.05_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")

# get summaries of FST distrs
summary(na.omit(maf05_all_to_all$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_cord_globmr$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_cord_globref$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_globmr_globref$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, all sample groups without WF03/1051
# import FST results
maf05_all_to_all_outl <- read.table("maf0.05_all_to_all_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf05_cord_globmr_outl <- read.table("maf0.05_cord_globmr_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf05_cord_globref_outl <- read.table("maf0.05_cord_globref_outl.weir.fst", header = TRUE, na.strings = "-nan")
maf05_globmr_globref_outl <- read.table("maf0.05_globmr_globref_outl.weir.fst", header = TRUE, na.strings = "-nan")

# get summaries of FST distrs
summary(na.omit(maf05_all_to_all_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_cord_globmr_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_cord_globref_outl$WEIR_AND_COCKERHAM_FST))
summary(na.omit(maf05_globmr_globref_outl$WEIR_AND_COCKERHAM_FST))
```

Fst distribution for all groups against each other:

**MAF = 0.00:**
| has WF03 | Group 1      | Group 2      | Min    | 1st Quart | Median | Mean   | 3rd Quart | Max  |
|----------|--------------|--------------|--------|-----------|--------|------- |-----------|------|
| yes      | all three    | N/A          | -0.106 | 0.0153    | 0.0829 | 0.170  | 0.248     | 1.00 |
| yes      | E. cordata   | E. glob intr | -0.676 | 0.0162    | 0.125  | 0.217  | 0.348     | 1.00 |
| yes      | E. cordata   | E. glob ref  | -1.00  | 0.0078    | 0.125  | 0.214  | 0.338     | 1.00 |
| yes      | E. glob intr | E. glob ref  | -0.195 | -0.0317   | -0.0133| 0.00966| 0.0385    | 0.833|
| no       | all three    | N/A          | -0.106 | 0.0153    | 0.0829 | 0.170  | 0.249     | 1.00 |
| no       | E. cordata   | E. glob intr | -0.676 | 0.0162    | 0.125  | 0.217  | 0.348     | 1.00 |
| no       | E. cordata   | E. glob ref  | -1.00  | 0.00784   | 0.125  | 0.214  | 0.339     | 1.00 |
| no       | E. glob intr | E. glob ref  | -0.196 | -0.0320   | -0.0136| 0.00966| 0.0391    | 0.916|

**MAF = 0.05:**
| has WF03 | Group 1      | Group 2      | Min    | 1st Quart | Median | Mean   | 3rd Quart | Max  |
|----------|--------------|--------------|--------|-----------|--------|------- |-----------|------|
| yes      | all three    | N/A          | -0.104 | 0.0422    | 0.155  | 0.234  | 0.355     | 1.00 |
| yes      | E. cordata   | E. glob intr | -0.676 | 0.0554    | 0.218  | 0.284  | 0.455     | 1.00 |
| yes      | E. cordata   | E. glob ref  | -0.392 | 0.0625    | 0.211  | 0.276  | 0.425     | 1.00 |
| yes      | E. glob intr | E. glob ref  | -0.137 | -0.0344   | -0.0138| 0.0103 | 0.0385    | 0.833|
| no       | all three    | N/A          | -0.106 | 0.0424    | 0.155  | 0.234  | 0.355     | 1.00 |
| no       | E. cordata   | E. glob intr | -0.676 | 0.0553    | 0.218  | 0.284  | 0.455     | 1.00 |
| no       | E. cordata   | E. glob ref  | -0.365 | 0.0625    | 0.211  | 0.276  | 0.425     | 1.00 |
| no       | E. glob intr | E. glob ref  | -0.178 | -0.0345   | -0.0140| 0.0103 | 0.0385    | 0.915|


**Conclusion:** Removing the strange samples doesn't seem to FST (and therefore relative variance in allele frequencies) between my sample groups.

## Pi and dXY

Also checked whether genome-wide sliding window distributions of nucleotide diversity (pi) and absolute divergence (dXY) were affected by the inclusion of WF03. To do so, calculated pi and dXY in sliding windows of 5000 bp with a 2500 bp step, excluding once sample from the _E. globulus_ reference group at a time.

```bash
# Performed on UFRC queue system. See stats_outlier_check.job for more details.
# Resources used: 400 Mb, 13 hrs

module load conda
ENV_DIR="/blue/soltis/kasey.pham/conda/envs"
conda activate euc_hyb_reseq

python stats_outlier_check.py
```

Then summarized distributions of stats for each jackknifed dataset in R.

e.g.,

```R
wdir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/outlier_check/pi_dxy"
infile_name <- "stats_outlier_check_no_WF03.tab"

setwd(wdir)
infile <- read.table(infile_name, header = TRUE, sep = "\t")

summary(infile$Pi)
summary(infile$Dxy)
```

**Results:**

**Nucleotide Diversity:**

| Sample Excluded | Accession | Min   | 1st Quart | Median | Mean  | 3rd Quart | Max    |
| --------------- | --------- | ----- | --------- | ------ | ----- | --------- | ------ |
| WF03            | 1051      | 0.00  | 4.51      | 15.73  | 20.26 | 30.19     | 199.15 |
| WA02            | 1119      | 0.00  | 4.53      | 15.77  | 20.30 | 30.24     | 199.20 |
| WB05            | 1081      | 0.00  | 4.45      | 15.43  | 19.98 | 29.72     | 198.76 |
| WC01            | 1129      | 0.00  | 4.27      | 14.79  | 19.30 | 28.62     | 192.42 |
| WC04            | 1043      | 0.00  | 4.53      | 15.76  | 20.29 | 30.24     | 198.83 |
| WD01            | 1130      | 0.00  | 4.53      | 15.77  | 20.29 | 30.22     | 199.01 |
| WD03            | 1116      | 0.00  | 4.52      | 15.75  | 20.28 | 30.22     | 198.73 |
| WF05            | 1059      | 0.00  | 4.52      | 15.76  | 20.28 | 30.23     | 197.54 |
| WH01            | 1056      | 0.00  | 4.52      | 15.76  | 20.28 | 30.22     | 194.44 |
| WH02            | 1128      | 0.00  | 4.41      | 15.40  | 19.86 | 29.60     | 194.22 |

**Absolute Divergence:**

| Sample Excluded | Accession | Min   | 1st Quart | Median | Mean | 3rd Quart | Max   |
| --------------- | --------- | ----- | --------- | ------ | ---- | --------- | ----- |
| WF03            | 1051      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WA02            | 1119      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.88  |
| WB05            | 1081      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WC01            | 1129      | 0.00  | 0.12      | 0.15   | 0.16 | 0.19      | 0.78  |
| WC04            | 1043      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WD01            | 1130      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WD03            | 1116      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.79  |
| WF05            | 1059      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WH01            | 1056      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |
| WH02            | 1128      | 0.00  | 0.13      | 0.16   | 0.17 | 0.20      | 0.78  |

**Conclusion:** There doesn't look to be an overall difference in at least the distribution for my diagnostic stats if I exclude WF03 versus excluding other samples. There could still be a spatial difference across the genome in how these stats are distributed... However, for now, I think it's acceptable to continue on including WF03 within my reference set since it's hard to tell where the differentiation is coming from.