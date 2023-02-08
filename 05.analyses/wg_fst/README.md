# Whole-genome Fst Calculations
Checked values of FST between populations with and without strange samples which stood out in whole genome PCA to try to get a better understanding of the distribution of variation in samples. WF03/1051 and WG04/2025 were distant in 2nd and 3rd PC plots from the rest of their sample groups, so did separate calculations of them against the different populations and ommitting them in different combinations.

```bash
module load vcftools/0.1.16

INFILE00="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.00/meehan_all_fil_maf0.00_snps.vcf"
INFILE05="/blue/soltis/kasey.pham/euc_hyb_reseq/call_snps/04.filter_snps/maf0.05/meehan_all_fil_maf0.05_snps.vcf"
LISTDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst"

# ALL SAMPLE GROUPS AGAINST EACH OTHER
# All, MAF=0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_all_to_all
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.00_cord_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_cord_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_globmr_globref

# All, MAF=0.05
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_all_to_all
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.05_cord_globmr
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_cord_globref
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.05_globmr_globref

# PROBLEMATIC SAMPLES AGAINST EACH SAMPLE GROUP
# WF03, MAF=0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_wf03_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.00_wf03_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.00_wf03_cord

# WG04, MAF=0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref.txt --out maf0.00_wg04_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --out maf0.00_wg04_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/WG04.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.00_wg04_cord

# WF03, MAF=0.05
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_wf03_globref
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR.txt --out maf0.05_wf03_globmr
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/WF03.txt --weir-fst-pop "$LISTDIR"/Ecordata.txt --out maf0.05_wf03_cord

# SAMPLE GROUPS AGAINST EACH OTHER WITH PROBLEMATIC SAMPLES REMOVED
# MAF = 0.00
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_all_to_all
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --out maf0.00_cord_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_cord_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_globmr_globref

# MAF = 0.05
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_all_to_all
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --out maf0.05_cord_globmr
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_cord_globref
vcftools --vcf "$INFILE05" --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.05_globmr_globref
```

Get genome-wide average in `R`:

```R
maf00_all_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.00_all"
maf05_all_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.05_all"
maf00_wf03_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.00_wf03"
maf00_wg04_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.00_wg04"
maf05_wf03_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.05_wf03"

# MAF=0.00, all sample groups against each other
setwd(maf00_all_dir)
maf00_all_to_all_tab <- read.table("maf0.00_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_all_to_all_tab$WEIR_AND_COCKERHAM_FST))
maf00_cord_globmr_tab <- read.table("maf0.00_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_cord_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_cord_globref_tab <- read.table("maf0.00_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_cord_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_globmr_globref_tab <- read.table("maf0.00_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_globmr_globref_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, all sample groups against each other
setwd(maf05_all_dir)
maf05_all_to_all_tab <- read.table("maf0.05_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_all_to_all_tab$WEIR_AND_COCKERHAM_FST))
maf05_cord_globmr_tab <- read.table("maf0.05_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_cord_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf05_cord_globref_tab <- read.table("maf0.05_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_cord_globref_tab$WEIR_AND_COCKERHAM_FST))
maf05_globmr_globref_tab <- read.table("maf0.05_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_globmr_globref_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.00, WF03 versus the three sample groups
setwd(maf00_wf03_dir)
maf00_wf03_globref_tab <- read.table("maf0.00_wf03_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_wf03_globmr_tab <- read.table("maf0.00_wf03_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_wf03_cord_tab <- read.table("maf0.00_wf03_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wf03_cord_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.00, WG04 versus the three sample groups
setwd(maf00_wg04_dir)
maf00_wg04_globref_tab <- read.table("maf0.00_wg04_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_wg04_globmr_tab <- read.table("maf0.00_wg04_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_wg04_cord_tab <- read.table("maf0.00_wg04_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_wg04_cord_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, WF03 versus the three sample groups
setwd(maf05_wf03_dir)
maf05_wf03_globref_tab <- read.table("maf0.05_wf03_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_globref_tab$WEIR_AND_COCKERHAM_FST))
maf05_wf03_globmr_tab <- read.table("maf0.05_wf03_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf05_wf03_cord_tab <- read.table("maf0.05_wf03_cord.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_wf03_cord_tab$WEIR_AND_COCKERHAM_FST))
```

Fst distribution for all groups against each other:
| MAF  | Group 1      | Group 2      | Min    | 1st Quart | Median | Mean   | 3rd Quart | Max  |
|------|--------------|--------------|--------|-----------|--------|------- |-----------|------|
| 0.00 | all three    | N/A          | -0.099 | 0.016     | 0.085  | 0.17   | 0.25      | 1.00 |
| 0.00 | E. cordata   | E. glob intr | -0.10  | 0.019     | 0.13   | 0.22   | 0.35      | 1.00 |
| 0.00 | E. cordata   | E. glob ref  | -0.15  | 0.0078    | 0.13   | 0.22   | 0.34      | 1.00 |
| 0.00 | E. glob intr | E. glob ref  | -0.11  | -0.031    | -0.013 | 0.0097 | 0.039     | 0.83 |
| 0.05 | all three    | N/A          | -0.099 | -0.043    | 0.16   | 0.23   | 0.36      | 1.00 |
| 0.05 | E. cordata   | E. glob intr | -0.10  | 0.055     | 0.22   | 0.28   | 0.018     | 1.00 |
| 0.05 | E. cordata   | E. glob ref  | -0.15  | 0.063     | 0.21   | 0.28   | 0.43      | 1.00 |
| 0.05 | E. glob intr | E. glob ref  | -0.14  | -0.034    | -0.014 | 0.010  | 0.039     | 0.83 |

Fst distribution for problematic samples against each sample group:

| MAF  | Group 1   | Group 2      | Min   | 1st Quart | Median | Mean  | 3rd Quart | Max  |
|------|-----------|--------------|-------|-----------|--------|------ |-----------|------|
| 0.00 | WF03/1051 | E. glob ref  | -1.09 | -0.36     | -0.22  | -0.11 | 0.056     | 1.00 |
| 0.00 | WF03/1051 | E. glob intr | -1.00 | -0.31     | -0.22  | -0.11 | 0.018     | 1.00 |
| 0.00 | WF03/1051 | E. cordata   | -1.09 | -0.20     | 0.051  | 0.19  | 0.72      | 1.00 |
| 0.00 | WG04/2025 | E. glob ref  | -1.09 | -0.36     | -0.26  | -0.13 | 0.00      | 1.00 |
| 0.00 | WG04/2025 | E. glob intr | -1.00 | -0.34     | -0.23  | -0.13 | -0.03     | 1.00 |
| 0.00 | WG04/2025 | E. cordata   | -1.19 | -0.22     | 0.052  | 0.19  | 0.76      | 1.00 |
| 0.05 | WF03/1051 | E. glob ref  | -1.09 | -0.36     | -0.18  | -0.094| 0.074     | 1.00 |
| 0.05 | WF03/1051 | E. glob intr | -1.00 | -0.28     | -0.16  | -0.077| 0.062     | 1.00 |
| 0.05 | WF03/1051 | E. cordata   | -1.09 | -0.099    | 0.17   | 0.27  | 0.82      | 1.00 |

Fst distribution for sample groups against each other with problematic samples omitted:


**Conclusion:** Divergence is pretty consistent between problem samples, their own group, and the other _E. globulus_ group but is much less than for problem samples with _E. cordata_. I can't see a reason to exclude them since they're exhibiting the pattern of divergence I'd expect of any other sample (can't speak to the degree of divergence though). HOWEVER: Since most of Fst values are negative (indicating larger variation within than between populations), it might not be the best metric to use for evaluating these samples. Will need to test whether results are robust to excluding WF03/1051 at the very least.