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
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_all_to_all
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --out maf0.00_cord_globmr
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Ecordata.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_cord_globref
vcftools --vcf "$INFILE00" --weir-fst-pop "$LISTDIR"/Eglobulus_MR_withoutWG04.txt --weir-fst-pop "$LISTDIR"/Eglobulus_ref_withoutWF03.txt --out maf0.00_globmr_globref

# MAF = 0.05
cd "$MAF05_DIR"
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

maf00_noWeird_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.00_noWeird"
maf05_noWeird_dir <- "/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/wg_fst/maf0.05_noWeird"

# ALL SAMPLE GROUPS 
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

# ALL SAMPLE GROUPS WITHOUT PROBLEMATIC SAMPLES
# MAF=0.00, all sample groups without WF03/1051 & WG04/2025
setwd(maf00_noWeird_dir)
maf00_all_to_all_tab <- read.table("maf0.00_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_all_to_all_tab$WEIR_AND_COCKERHAM_FST))
maf00_cord_globmr_tab <- read.table("maf0.00_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_cord_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf00_cord_globref_tab <- read.table("maf0.00_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_cord_globref_tab$WEIR_AND_COCKERHAM_FST))
maf00_globmr_globref_tab <- read.table("maf0.00_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf00_globmr_globref_tab$WEIR_AND_COCKERHAM_FST))

# MAF=0.05, all sample groups without WF03/1051 & WG04/2025
setwd(maf05_noWeird_dir)
maf05_all_to_all_tab <- read.table("maf0.05_all_to_all.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_all_to_all_tab$WEIR_AND_COCKERHAM_FST))
maf05_cord_globmr_tab <- read.table("maf0.05_cord_globmr.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_cord_globmr_tab$WEIR_AND_COCKERHAM_FST))
maf05_cord_globref_tab <- read.table("maf0.05_cord_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_cord_globref_tab$WEIR_AND_COCKERHAM_FST))
maf05_globmr_globref_tab <- read.table("maf0.05_globmr_globref.weir.fst", header = TRUE, na.strings = "-nan")
summary(na.omit(maf05_globmr_globref_tab$WEIR_AND_COCKERHAM_FST))
```

Fst distribution for all groups against each other:
| MAF  | Group 1      | Group 2      | Min    | 1st Quart | Median | Mean   | 3rd Quart | Max  |
|------|--------------|--------------|--------|-----------|--------|------- |-----------|------|
| 0.00 | all three    | N/A          | -0.099 | 0.016     | 0.085  | 0.17   | 0.25      | 1.00 |
| 0.00 | E. cordata   | E. glob intr | -0.10  | 0.019     | 0.13   | 0.22   | 0.35      | 1.00 |
| 0.00 | E. cordata   | E. glob ref  | -0.15  | 0.0078    | 0.13   | 0.22   | 0.34      | 1.00 |
| 0.00 | E. glob intr | E. glob ref  | -0.11  | -0.031    | -0.013 | 0.0097 | 0.039     | 0.83 |
| 0.05 | all three    | N/A          | -0.099 | -0.043    | 0.16   | 0.23   | 0.36      | 1.00 |
| 0.05 | E. cordata   | E. glob intr | -0.10  | 0.055     | 0.22   | 0.28   | 0.45      | 1.00 |
| 0.05 | E. cordata   | E. glob ref  | -0.15  | 0.063     | 0.21   | 0.28   | 0.43      | 1.00 |
| 0.05 | E. glob intr | E. glob ref  | -0.14  | -0.034    | -0.014 | 0.010  | 0.039     | 0.83 |

Fst distribution for sample groups against each other with WF03/1051 and WG04/2025 omitted:

| MAF  | Group 1      | Group 2      | Min    | 1st Quart | Median | Mean   | 3rd Quart | Max  |
|------|--------------|--------------|--------|-----------|--------|------- |-----------|------|
| 0.00 | all three    | N/A          | -0.10  | 0.016     | 0.086  | 0.17   | 0.25      | 1.00 |
| 0.00 | E. cordata   | E. glob intr | -0.096 | 0.020     | 0.13   | 0.22   | 0.35      | 1.00 |
| 0.00 | E. cordata   | E. glob ref  | -0.15  | 0.0078    | 0.13   | 0.22   | 0.34      | 1.00 |
| 0.00 | E. glob intr | E. glob ref  | -0.11  | -0.032    | -0.014 | 0.0098 | 0.040     | 0.91 |
| 0.05 | all three    | N/A          | -0.10  | -0.044    | 0.16   | 0.24   | 0.35      | 1.00 |
| 0.05 | E. cordata   | E. glob intr | -0.096 | 0.055     | 0.22   | 0.28   | 0.45      | 1.00 |
| 0.05 | E. cordata   | E. glob ref  | -0.15  | 0.066     | 0.21   | 0.28   | 0.43      | 1.00 |
| 0.05 | E. glob intr | E. glob ref  | -0.11  | -0.035    | -0.014 | 0.010  | 0.039     | 0.91 |


**Conclusion:** Removing the strange samples doesn't seem to affect differentiation measures (and therefore relative variance in allele frequencies) between my sample groups.

## Nei's D
