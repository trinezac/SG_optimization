library(vegan)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

load("/Users/trinezachariasen/Downloads/bakhed/abundance_phyloseq.RData")

meta <- read.table("/Users/trinezachariasen/Downloads/bakhed/bakhed_reads_meta.txt", header = TRUE, sep = "\t")

# removing the samples with no metadata
no_meta_samples <-c("ERR525840","ERR525841","ERR525842","ERR525843","ERR525976","ERR525977","ERR525978","ERR525979")
physeq_pruned <- prune_samples(!sample_names(final.physeq) %in% no_meta_samples, final.physeq)

# create beta diversity matrix
beta_mat <- distance(physeq_pruned, method = "bray")

# perform classical multidimensional scaling
pcoa <- cmdscale(beta_mat, eig=TRUE)

eig <- pcoa$eig
# Calculate total sum of eigenvalues
total_eig <- sum(eig)
# Calculate percentage explained by each PC
pc1_perc_expl <- round(eig[1] / total_eig * 100, 2)
pc2_perc_expl <- round(eig[2] / total_eig * 100, 2)

# convert to data frame for plotting
pcoa_df <- data.frame(pcoa$points, row.names=rownames(pcoa$points))
colnames(pcoa_df) <- c("PC1", "PC2")

# Adding the stage/sample type to the dataframe
stage <- c()
for (entry in rownames(pcoa_df)){
  entry_stage <- meta$stage[meta$sample==entry]
  stage <- c(stage, entry_stage)
}
pcoa_meta <- cbind(pcoa_df, stage)

# plot PCoA using ggplot2
pcoa_plot <- ggplot(pcoa_meta, aes(x=PC1, y=PC2, color=stage)) + 
  geom_point(size=1) + stat_ellipse() +
  labs(x=paste("PC1 ", pc1_perc_expl, "%", sep=""), y=paste("PC2 ", pc2_perc_expl, "%", sep=""), color="Stage")
pcoa_plot


permanova <- adonis(beta_mat ~ stage, method="bray", permutations=999)
# P-value
permanova$aov.tab

## Regenerating their heatmap

# Table S5, sheet 2: Signature genera at each stage for vaginally delivered infants
vag_taxa_list <- c("Enterococcus", "Escherichia", "Escherichia/Shigella", "Gemella", "Propionibacterium", "Rothia", "Staphylococcus", "Streptococcus", "Actinomyces", "Anaerococcus", "Bifidobacterium", "Collinsella", "Eggerthella", "Granulicatella", "Lactobacillus", "Mobiluncus", "Peptoniphilus", "Veillonella", "Aggregatibacter", "Anaerostipes", "Anaerotruncus", "Bacteroides", "Clostridiales", "Clostridium", "Eikenella", "Fusobacterium", "Lachnospiraceae", "Megasphaera", "Roseburia", "Ruminococcus", "Acidaminococcus", "Akkermansia", "Alistipes", "Bilophila", "Blautia", "Bryantella", "Butyrivibrio", "Capnocytophaga/Paraprevotella", "Coprobacillus", "Coprococcus", "Desulfovibrio", "Dialister", "Dorea", "Erysipelotrichaceae", "Eubacterium", "Faecalibacterium", "Holdemania", "Leuconostoc", "Methanobrevibacter", "Odoribacter", "Parasutterella", "Phascolarctobacterium", "Porphyromonas", "Prevotella", "Ruminococcaceae", "Subdoligranulum", "Turicibacter", "Victivallis")

# Table S5, sheet 4: Signature genera at each stage for C-section delivered infants.
csec_taxa_list <- c("Aggregatibacter", "Gemella", "Haemophilus", "Rothia", "Staphylococcus","Actinomyces", "Bifidobacterium", "Anaerostipes", "Anaerotruncus","Blautia", "Clostridiales", "Coprobacillus", "Eggerthella", "Eubacterium","Lachnospiraceae", "Prevotella", "Roseburia", "Ruminococcus", "Turicibacter","Alistipes", "Bacteroides", "Bilophila", "Bryantella", "Butyrivibrio","Coprococcus", "Dialister", "Dorea", "Erysipelotrichaceae", "Faecalibacterium","Holdemania", "Leuconostoc", "Methanobrevibacter", "Parasutterella","Ruminococcaceae", "Subdoligranulum", "Sutterella", "Victivallis")


tax_df <- as.data.frame(tax_table(final.physeq))
abundance_df <- as.data.frame(otu_table(final.physeq))

species_df <- abundance_df
species <- c()
for (cluster in colnames(abundance_df)){
  cluster_species <- tax_df[cluster,"Genus"]
  species <- c(species, cluster_species)
}

colnames(species_df) <- species

# Aggregate columns by name and sum them up
species_summed <- sapply(split.default(species_df, names(species_df)), rowSums, na.rm = TRUE)

t_species <- t(species_summed)
t_species <- t_species[, !(colnames(t_species) %in% no_meta_samples)] # removing the samples with no metadata

vag_tax_present <- rownames(t_species)[rownames(t_species)%in% vag_taxa_list]
ordered_vag_tax <- rev(vag_taxa_list[vag_taxa_list %in% vag_tax_present])

t_species<-t_species[rownames(t_species)%in% vag_taxa_list,]
t_species <- t_species[match(ordered_vag_tax, rownames(t_species)), ] # ordering according to list in TS5

t_species <- t_species[, (colnames(t_species) %in% meta$sample[meta$Delivery.mode=="vaginal"])]

#deselecting the sample_ids, which did not assemble and also the samples with no metadata
no_meta_samples <-c("ERR525840","ERR525841","ERR525842","ERR525843","ERR525976","ERR525977","ERR525978","ERR525979")
pruned_meta <- meta[!meta$sample %in% no_meta_samples,]
samples_new <- pruned_meta$sample[pruned_meta$stage=="new"][pruned_meta$sample[pruned_meta$stage=="new"] %in% colnames(t_species)]
samples_4m <- pruned_meta$sample[pruned_meta$stage=="4m"][pruned_meta$sample[pruned_meta$stage=="4m"] %in% colnames(t_species)]
samples_12m <- pruned_meta$sample[pruned_meta$stage=="12m"][pruned_meta$sample[pruned_meta$stage=="12m"] %in% colnames(t_species)]
samples_mom <- pruned_meta$sample[pruned_meta$stage=="Mom"][pruned_meta$sample[pruned_meta$stage=="Mom"] %in% colnames(t_species)]

# generating abund df for each stage
abund_new <- t_species[,samples_new]
abund_4m <- t_species[,samples_4m]
abund_12m <- t_species[,samples_12m]
abund_mom <- t_species[,samples_mom]

sorted_abund <- cbind(abund_new, abund_4m, abund_12m, abund_mom)
n_samples <- c(length(samples_new), length(samples_4m),length(samples_12m),length(samples_mom))
mydf <- data.frame(row.names =colnames(sorted_abund), Stage = c(rep("New", n_samples[1]), rep("4m",n_samples[2]), rep("12m", n_samples[3]), rep("Mom",n_samples[4])))

pheatmap(log(sorted_abund+0.01), cluster_cols = F, cluster_rows = F, annotation_col = mydf, show_rownames = T, show_colnames = F, gaps_col = cumsum(n_samples), gaps_row=c(19,25,32))

#pheatmap((sorted_abund), cluster_cols = F, cluster_rows = F, annotation_col = mydf, show_rownames = T, show_colnames = F, gaps_col = cumsum(n_samples), gaps_row=c(19,25,32))

####### For csection

t_species <- t(species_summed)
t_species <- t_species[, !(colnames(t_species) %in% no_meta_samples)] # removing the samples with no metadata

csec_tax_present <- rownames(t_species)[rownames(t_species)%in% csec_taxa_list]
ordered_vag_tax <- rev(csec_taxa_list[csec_taxa_list %in% csec_tax_present])

t_species<-t_species[rownames(t_species)%in% csec_taxa_list,]
t_species <- t_species[match(ordered_vag_tax, rownames(t_species)), ] # ordering according to list in TS5

t_species <- t_species[, (colnames(t_species) %in% meta$sample[meta$Delivery.mode=="section"])]

#deselecting the sample_ids, which did not assemble and also the samples with no metadata
no_meta_samples <-c("ERR525840","ERR525841","ERR525842","ERR525843","ERR525976","ERR525977","ERR525978","ERR525979")
pruned_meta <- meta[!meta$sample %in% no_meta_samples,]
samples_new <- pruned_meta$sample[pruned_meta$stage=="new"][pruned_meta$sample[pruned_meta$stage=="new"] %in% colnames(t_species)]
samples_4m <- pruned_meta$sample[pruned_meta$stage=="4m"][pruned_meta$sample[pruned_meta$stage=="4m"] %in% colnames(t_species)]
samples_12m <- pruned_meta$sample[pruned_meta$stage=="12m"][pruned_meta$sample[pruned_meta$stage=="12m"] %in% colnames(t_species)]
samples_mom <- pruned_meta$sample[pruned_meta$stage=="Mom"][pruned_meta$sample[pruned_meta$stage=="Mom"] %in% colnames(t_species)]

# generating abund df for each stage
abund_new <- t_species[,samples_new]
abund_4m <- t_species[,samples_4m]
abund_12m <- t_species[,samples_12m]
abund_mom <- t_species[,samples_mom]

sorted_abund <- cbind(abund_new, abund_4m, abund_12m, abund_mom)
n_samples <- c(length(samples_new), length(samples_4m),length(samples_12m),length(samples_mom))
mydf <- data.frame(row.names =colnames(sorted_abund), Stage = c(rep("New", n_samples[1]), rep("4m",n_samples[2]), rep("12m", n_samples[3]), rep("Mom",n_samples[4])))

pheatmap(log(sorted_abund+0.01), cluster_cols = F, cluster_rows = F, annotation_col = mydf, show_rownames = T, show_colnames = F, gaps_col = cumsum(n_samples), gaps_row=c(10,19,21),border_color = NA)






