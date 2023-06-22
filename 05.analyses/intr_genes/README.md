# Genes in Genome Outlier Windows

## Get genes within pi/dxy outlier windows

Got genes which fall within outlier windows using [`BEDTools`](https://bedtools.readthedocs.io).

```bash
module load bedtools/2.30.0

BED_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/genome_scan"
ANNOT_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46"
declare -a NAMELIST=(pi_dxy_outl_p05 pi_dxy_outl_p10 pi_dxy_outl_p15)

for NAME in "${NAMELIST[@]}"
do
    bedtools intersect -a "$BED_DIR"/"$NAME".bed -b "$ANNOT_DIR"/EGLOB-X46.v1.0.annotation.gff -wb > "$NAME"_genes.bed
done
```

Retrieved names of genes in overlapping regions.

```bash
declare -a NAMELIST=(pi_dxy_outl_p05 pi_dxy_outl_p10 pi_dxy_outl_p15)
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"
module load python/3.8

for NAME in "${NAMELIST[@]}"
do
    python "$SCRIPTS_DIR"/parse_intersects.py "$NAME"_genes.bed "$NAME"_genelist.txt
done
```

Extracted genes from reference matching genes in overlap regions.

```R
library(seqinr)

ref_genes <- read.fasta("/blue/soltis/kasey.pham/euc_hyb_reseq/refs/Eglobulus_genome_X46/EGLOB-X46.v1.0.genes.fa")
p05_list <- read.table("pi_dxy_outl_p05_genelist.txt", header = FALSE, as.is = TRUE)$V1
p10_list <- read.table("pi_dxy_outl_p10_genelist.txt", header = FALSE, as.is = TRUE)$V1
p15_list <- read.table("pi_dxy_outl_p15_genelist.txt", header = FALSE, as.is = TRUE)$V1

# remove duplicate genes from cutoff threshold group before
p15_list <- p15_list[which(!p15_list %in% p10_list)]
p10_list <- p10_list[which(!p10_list %in% p05_list)]

p05_genes <- ref_genes[which(names(ref_genes) %in% p05_list)]
p10_genes <- ref_genes[which(names(ref_genes) %in% p10_list)]
p15_genes <- ref_genes[which(names(ref_genes) %in% p15_list)]

write.fasta(p05_genes, names(p05_genes), "pi_dxy_outl_p05_genes.fasta")
write.fasta(p10_genes, names(p10_genes), "pi_dxy_outl_p10_genes.fasta")
write.fasta(p15_genes, names(p15_genes), "pi_dxy_outl_p15_genes.fasta")
```

## Get putative function

Examined annotation notes from BED files of genes in overlap windows.

**5%:**
| Chr   | Start    | End      | X46 Gene ID | Annotation                                                                              |
| ----- | -------- | -------- | ----------- |---------------------------------------------------------------------------------------- |
| Chr08 | 67870038 | 67870334 | ANN14063    | Similar to PDF2.5: Defensin-like protein 6 (Arabidopsis thaliana OX=3702)               |
| Chr08 | 67872593 | 67872888 | ANN14064    | Similar to PDF2.5: Defensin-like protein 6 (Arabidopsis thaliana OX=3702)               |

Seem too short to be real genes; Look into further.

**10%:**

| Chr   | Start    | End      | X46 Gene ID | Annotation                                                                              |
| ----- | -------- | -------- | ----------- |---------------------------------------------------------------------------------------- |
| Chr02 | 14118055 | 14119080 | ANN17914    | Similar to BG5: Probable glucan endo-1,3-beta-glucosidase BG5 (A. thaliana OX=3702)     |
| Chr03 | 17192399 | 17195067 | ANN24059    | Protein of unknown function                                                             |
| Chr03 | 17195990 | 17197937 | ANN24060    | Protein of unknown function                                                             |
| Chr03 | 57876689 | 57883889 | ANN22523    | Similar to MPP: Mitochondrial-processing peptidase subunit alpha (S. tuberosum OX=4113) |
| Chr06 | 16121099 | 16123567 | ANN15393    | Similar to BRN1: Protein BEARSKIN1 (Arabidopsis thaliana OX=3702)                       |

**15%:**
| Chr   | Start    | End      | X46 Gene ID | Annotation                                                                              |
| ----- | -------- | -------- | ----------- |---------------------------------------------------------------------------------------- |
| Chr03 | 29220495 | 29235591 | ANN21306    | Similar to SPCC1672.07: U3 small nucleolar RNA-associated protein 21 homolog (S. pombe (strain 972 / ATCC 24843) OX=284812) |
| Chr06 | 4009547  | 4015598  | ANN14454    | Similar to PEPKR2: Serine/threonine-protein kinase PEPKR2 (A. thaliana OX=3702)         |
| Chr06 | 17411506 | 17417049 | ANN15492    | Similar to PSI2: Protein PSK SIMULATOR 2 (Arabidopsis thaliana OX=3702)                 |
| Chr06 | 23955855 | 23960369 | ANN15818    | Similar to GAMMACA1: Gamma carbonic anhydrase 1, mitochondrial (A. thaliana OX=3702)    |
| Chr06 | 46336952 | 46342707 | ANN17317    | Similar to OSCBPY: Beta-amyrin synthase (Betula platyphylla OX=78630)                   |
| Chr06 | 46333564 | 46335504 | ANN17316    | Protein of unknown function                                                             |
| Chr08 | 6198010  | 6202241  | ANN10679    | Similar to SEOB: Protein SIEVE ELEMENT OCCLUSION B (Arabidopsis thaliana OX=3702)       |
| Chr09 | 36559718 | 36561710 | ANN30585    | Similar to UGT85A24: 7-deoxyloganetin glucosyltransferase (Gardenia jasminoides OX=114476) |