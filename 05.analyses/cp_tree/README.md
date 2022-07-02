# Chloroplast Haplotype Analysis

## Get chloroplast alignment

**Create FASTA files with all chloroplast region assemblies and reformat:**
Wrote python scripts to concatenate chloroplast assemblies from each sample into separate FASTA files for each region (files currently hardcoded!) and reformat FASTA files to split sequences into multiple lines.

```bash
module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPTS_DIR"/make_cp_region_fas.py
# split each sequence into multiple lines
python "$SCRIPTS_DIR"/reformat_fasta.py 80 irb_temp.fas irb_assemblies.fas
python "$SCRIPTS_DIR"/reformat_fasta.py 80 lsc_temp.fas lsc_assemblies.fas
python "$SCRIPTS_DIR"/reformat_fasta.py 80 ssc_temp.fas ssc_assemblies.fas
rm *_temp.fas
```

**Align chloroplast regions:**
```bash
# Run via job queue on UFRC, see mafft_ssc.job, mafft_irb.job, mafft_lsc.job for more details.
# Resources used: ~500 Mb, 11 min maximum.
# One example of the three regions shown below.

module load mafft/7.490

mafft --auto --thread 16 lsc_assemblies.fas > lsc_assemblies_aligned.fas
```

**Manual editing of alignments:**
* Manually adjusted a few bases among several samples in lsc alignment in `BioEdit`. Saved as `lsc_assemblies_adjusted.fas`
* Manually adjusted a few bases among several samples in ssc alignment in `BioEdit`. Saved as `ssc_assemblies_adjusted.fas`


**Concatenate regions:**
```bash
module load python/3.8
SCRIPTS_DIR="/blue/soltis/kasey.pham/euc_hyb_reseq/scripts"

python "$SCRIPTS_DIR"/concatenate_fasta.py irb_assemblies_aligned.fas,lsc_assemblies_adjusted.fas,ssc_assemblies_adjusted.fas concatenated_cp_aligned.fas
```

## Chloroplast tree

**Phylogenetic analysis in [`IQTree`](http://www.iqtree.org/):**
```bash
# Run via job queueing in UFRC; see iqtree.job for more details
# Resources: 2 Mb, 15 sec

module load iq-tree/2.1.3
WDIR="/blue/soltis/kasey.pham/euc_hyb_reseq/analyses/cp_tree"

iqtree -s "$WDIR"/concatenated_cp_aligned.fas -st DNA -B 1000 -m MFP -nt 12 -pre concatenated_cp_aligned
```

**Plot phylogeny in `R`:**

```R
# Done on local computer
library(ape)
# init file names
meta_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
intree_name <- "concatenated_cp_aligned.treefile"

# read in files
meta_table <- read.csv(paste(meta_name, sep = ""), header = TRUE, as.is = TRUE)
intree <- read.tree(intree_name)

# add outgroup labels
meta_table <- rbind(meta_table, c("HM347959.1", "", "E. grandis"))
meta_table <- rbind(meta_table, c("KC180790.1", "", "E. saligna"))

# replace tip names
label_order <- match(intree$tip.label, meta_table$RAPiD_ID)
replacemt_labels <- paste(meta_table[label_order, "Accession"], meta_table[label_order, "Taxon"], sep = "_")
intree$tip.label <- replacemt_labels

# root tree
intree <- root(intree, c("E. grandis_"), resolve.root = TRUE)

# plot tree
tip_colors <- sapply(meta_table[label_order, "Taxon"], function(x) ifelse(x == "cord_MR", "goldenrod1", ifelse(x == "glob_MR", "black", "deepskyblue4")))

plot(intree, tip.color = tip_colors, cex = 1, show.node.label = TRUE, adj = 1, label.offset = 0.5)
```

![chloroplast phylogeny results, tips labeled by accession and sample species](cp_tree_rough.png "Chloroplast Phylogeny")