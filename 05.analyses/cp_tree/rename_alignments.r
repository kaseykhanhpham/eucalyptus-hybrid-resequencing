library(seqinr)
# init file names
meta_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/03.seq_analysis/sample_spp_table.csv"
cp_hap_meta_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/00.metadata/01.field_sampling/40samples_KaseyUF_updateDec2020.csv"
align_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/cp_tree/alignments/irb_assemblies_aligned.fas"
output_name <- "C:/Users/Kasey/OneDrive - University of Florida/Grad School Documents/Projects/eucalyptus-hybrid-resequencing/05.analyses/cp_tree/alignments/irb_assemblies_aligned_renamed.fas"

# Read inputs
meta_table <- read.csv(meta_name, header = TRUE, as.is = TRUE)
cp_table <- read.csv(cp_hap_meta_name, header = TRUE, as.is = TRUE)
alignmt <- read.fasta(align_name)

# Add outgroups
meta_table <- rbind(meta_table, c("HM347959.1", "", "E. grandis"))
meta_table <- rbind(meta_table, c("KC180790.1", "", "E. saligna"))

# make new labels
label_order <- match(names(alignmt), meta_table$RAPiD_ID)
cp_order <- match(meta_table[label_order, "Accession"], cp_table$RJgeno)
replacemt_labels <- paste(meta_table[label_order, "Accession"], meta_table[label_order, "Taxon"], cp_table[cp_order, "JLA."], sep = "_")

# re-write names and write to file
names(alignmt) <- replacemt_labels
write.fasta(alignmt, names = replacemt_labels, output_name)