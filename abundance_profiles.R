# Converting the read count matrix to abundance profiles

library(phyloseq)
library(forcats)
library(ggplot2)
library(MASS)
library(ggpubr) # for plotting the rel abundance together
library(DAtest)
require(dplyr)
require(tibble)
library(gridExtra)

work_dir = "/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/Github_versions/"
GeneLengths <- readRDS(paste(work_dir, 'SGC_genelengths.RDS', sep="")) # The gene lengths
outputs_parallel <- readRDS(file=paste(work_dir,'MSP_SGC_refined_coabundance.RDS',sep='')) # contain the SGs of the Clusters
Clusterlist <- readRDS(paste(work_dir, "MSP_coabundance_SGC.RDS", sep="")) # read count matrices for the Clusters
tax_df <- readRDS(paste(work_dir,"tax_df.RDS", sep=""))
true_tax <- t(read.csv("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/true_abundance.csv", header = FALSE, sep=";", row.names=1, col.names = c("",rownames(init.read.matrix))))


#setting important variables
gene_index <- seq(1,length(GeneLengths))
gene_names <- names(GeneLengths)
n.mapped.minimum <- 3 #The number of genes that needs reads that map to count the cluster as present
n.genes <- 100 # number of signature genes

# inserting NA for the Clusters and levels that do not have an annotation
taxmat <- matrix("NA", nrow = length(names(Clusterlist)), ncol = 8)
taxmat <- matrix(0, nrow = length(names(Clusterlist)), ncol = 8)

colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
rownames(taxmat) <- names(Clusterlist)

for (cluster in names(Clusterlist)){
  tax <- tax_df$Taxonomy[tax_df$MSP==cluster]
  tax <- strsplit(tax, "strain")[[1]]
  row <- rep(0,8)
  row[7:(6+length(tax))] <- as.character(tax)
  taxmat[cluster,] <- row
}

# The readcounts are normalized according to gene lengths
init.read.matrix <- matrix(NA, nrow=dim(Clusterlist[[1]])[2], ncol=length(Clusterlist))
final.read.matrix <- matrix(NA, nrow=dim(Clusterlist[[1]])[2], ncol=length(Clusterlist))
sample.names <- colnames(Clusterlist[[1]])
rownames(init.read.matrix) <- sample.names # setting the sample names as rownames
rownames(final.read.matrix) <- sample.names
colnames(init.read.matrix) <- names(Clusterlist) # setting the cluster ID's as colnames
colnames(final.read.matrix) <- names(Clusterlist)
init.Clusterlist <- Clusterlist
final.Clusterlist <- Clusterlist


for (id in names(Clusterlist)){
  #removing the samples, with less than 3 genes with mapped reads
  init.gene.names <- outputs_parallel[[id]]$genes$init
  init.colsum <- colSums(Clusterlist[[id]][init.gene.names, ])
  init.mapped <- names(init.colsum[init.colsum >= n.mapped.minimum])
  init.not_mapped <- sample.names[!sample.names %in% init.mapped]
  init.Clusterlist[[id]][,init.not_mapped] <- 0 # setting the counts to 0 if less than n.mapped.minimum genes have reads that map

  # The readcounts are divided by the gene length
  init.reads <- init.Clusterlist[[id]][init.gene.names, ] / GeneLengths[init.gene.names]

  # summing the read counts for the id/cluster/MGS
  init.read.matrix[, id] <- colSums(init.reads)

  # # Repeat for the final SG
  final.gene.names <- outputs_parallel[[id]]$genes$best
  final.colsum <- colSums(Clusterlist[[id]][final.gene.names, ])
  final.mapped <- names(final.colsum[final.colsum >= n.mapped.minimum])
  final.not_mapped <- sample.names[!sample.names %in% final.mapped]
  final.Clusterlist[[id]][,final.not_mapped] <- 0 # setting the counts to 0 if less than n.mapped.minimum genes have reads that map

  # The readcounts are divided by the gene length
  final.reads <- final.Clusterlist[[id]][final.gene.names, ] / GeneLengths[final.gene.names]

  # summing the read counts for the id/cluster/MGS
  final.read.matrix[, id] <- colSums(final.reads)
}

init.abundance <- init.read.matrix/rowSums(init.read.matrix)
final.abundance <- final.read.matrix/rowSums(final.read.matrix)


####################### Visualizing #######################

init.otu.table <- otu_table(init.abundance, taxa_are_rows=FALSE) 
final.otu.table <- otu_table(final.abundance, taxa_are_rows = FALSE)
tax.table <- tax_table(taxmat)

init.physeq <- phyloseq(init.otu.table, tax.table)
final.physeq <-  phyloseq(final.otu.table, tax.table)

#save(init.physeq, file = paste(data_dir, "SGC_phyloseq_init.RData", sep = ""))
#save(final.physeq, file = paste(data_dir, "SGC_phyloseq_refined.RData", sep = ""))

#load(paste(data_dir, "SGC_phyloseq_init.RData", sep = ""))
#load(paste(data_dir, "SGC_phyloseq_refined.RData", sep = ""))
col_names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")

#sorting according to phyla
i.ps_species <- tax_glom(init.physeq, taxrank = c("Species"))#, "Class", NArm = TRUE)
taxa_names(i.ps_species) <- tax_table(i.ps_species)[,"Species"]
ps_species <- tax_glom(final.physeq, taxrank = c("Species"))
taxa_names(ps_species) <- tax_table(ps_species)[, "Species"]

init.ps_species <- psmelt(i.ps_species)
final.ps_species <- psmelt(ps_species)

# displaying the abundance of each class in a separate figures
# psmelt(i.ps_species) %>%
#   ggplot(data = ., aes(x = Sample, y = Abundance)) +
#   geom_boxplot(outlier.shape  = NA) +
#   geom_jitter(aes(color = OTU), height = 0, width = .2) +
#   labs(x = "", y = "Abundance\n", title = "Sample-specific relative abundance of each taxonomic species of the initial SGs from MSPminer") +
#   facet_wrap(~ OTU, scales = "free") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none")
# 
# 
# psmelt(ps_species) %>%
#   ggplot(data = ., aes(x = Sample, y = Abundance)) +
#   geom_boxplot(outlier.shape  = NA) +
#   geom_jitter(aes(color = OTU), height = 0, width = .2) +
#   labs(x = "Samples", y = "Abundance\n", title = "Sample-specific relative abundance of each taxonomic species of the Refined SG") +
#   facet_wrap(~ OTU, scales = "free") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.title = element_blank(),legend.position="none")
# 

# Displaying the relative abundance as barplot
col <- c("#9ebcda",'#de77ae','#33a02c','#a6cee3','#a6cee3','#8dd3c7', '#e31a1c','#fdbf6f','#dd1c77','#ff7f00','#b3de69','#fb9a99','#cab2d6','#6a3d9a',"#9ebcda",'#de77ae','#33a02c','#a6cee3')

vulgatus <- psmelt(final.physeq)[psmelt(final.physeq)$Species== "Bacteroides vulgatus",]
lvls <- vulgatus[order(vulgatus$Abundance, decreasing = TRUE), c(2)]

filling <- fct_reorder(psmelt(final.physeq)$Species, psmelt(final.physeq)$Abundance) #reuse the colors for both figures

p1 <- psmelt(final.physeq) %>% 
  ggplot(data=.,  aes(x=substring(Sample, 8), y = Abundance, fill = filling)) + 
  geom_bar(stat = "identity") +
  xlab("Sample") + theme_minimal() +
  ylab("Species-level relative abundance (%)") + 
  ggtitle("B: Refined Signature Genes") + #, subtitle = "Each bar represents a sample. Samples are arranged by increasing abundance of 'Ruminococcus sp. SR1/5'.")+
#  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size = 11)) +
  labs(fill='Species') + scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(0, 1), expand = c(0, 0)) 

p2 <- psmelt(init.physeq) %>%
  ggplot(data=.,  aes(x=substring(Sample, 8), y = Abundance, fill = filling)) + 
  geom_bar(stat = "identity") +
  xlab("Sample") +
  ylab("Species-level relative abundance (%)") + 
  ggtitle("A: Signature Genes from MSPminer") + #, subtitle = "Each bar represents a sample. Samples are arranged by increasing abundance of 'Ruminococcus sp. SR1/5'.")+
#  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size = 11)) +
  labs(fill='Species') + scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(0, 1), expand = c(0, 0))


abundance_plots <- ggarrange(p2, p1 + theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        axis.title.y = element_blank() ), ncol=2, nrow=1, widths=c(1,0.91), common.legend=TRUE, legend = "bottom")#  labels=c("A", "B"), 

annotate_figure(abundance_plots, top = text_grob("Sample-specific relative abundance of species-level taxa", 
                                      color = "black", face = "bold", size = 15))


##############################################################################################################################################
# We want to compare the initial and the refined relative abundances with the TRUE abundance
##############################################################################################################################################


# The MSP's with identical taxa are collapsed into a single entity
init_tax_abund <- matrix(0, nrow = dim(true_tax)[1], ncol = dim(true_tax)[2], dimnames=list(row.names(true_tax), colnames(true_tax)))
final_tax_abund <- matrix(0, nrow = dim(true_tax)[1], ncol = dim(true_tax)[2], dimnames=list(row.names(true_tax), colnames(true_tax)))
for (j in 1:length(colnames(true_tax))){
  col = colnames(true_tax)[j]
  col = stringr::str_replace_all(col, "\\s", " ")  #changing the encoding
  for (i in 1:length(tax_df$Taxonomy)){
    if (col == tax_df$Taxonomy[i]){
      init_tax_abund[,j] <- init_tax_abund[j] + init.read.matrix[,tax_df$MSP[i]]
      final_tax_abund[,j] <- final_tax_abund[j] + final.read.matrix[,tax_df$MSP[i]]
    }
  }
}

# Normalizing
init_tax_abund_norm <- ((init_tax_abund)/rowSums(init_tax_abund)*100)
final_tax_abund_norm <- ((final_tax_abund)/rowSums(final_tax_abund)*100)

# Subtracting the true abundance from the predicted 
init_true <- init_tax_abund_norm-true_tax
final_true <- final_tax_abund_norm-true_tax

# Testing whether the predicted relative abundance distributions are identical
wilcox.test(init_true, final_true, paired = TRUE, conf.int=TRUE)
wilcox_effsize(data.frame("init_true"=as.vector(init_true), "final_true"=as.vector(final_true)), init_true ~ final_true) # takes a long time to run

# Visualizing
final_true_sum <- colSums(final_true)
init_true_sum <- colSums(init_true)

violin.plot.dat <- as.data.frame(rbind(cbind(colSums(init_true),rep('MSPminer',length(colSums(init_true)))),cbind(colSums(final_true),rep('Refined',length(colSums(final_true))))))
p <- ggplot(violin.plot.dat, aes(x=V2, y=as.numeric(V1), fill=V2)) + 
  geom_violin(trim=FALSE)  + theme_minimal() + geom_boxplot(width=0.1, fill="white") + scale_fill_brewer(palette="Blues") +
  labs(title="Plot of divergence between the true and predicted abundance", x="Method", y = "Difference in abundance (method-true)",fill="Method") 
p

# library(beeswarm)
# 
# beeswarm(abs(as.numeric(violin.plot.dat$V1[violin.plot.dat$V2 == "MSPminer"])) ~ abs(as.numeric(violin.plot.dat$V1[violin.plot.dat$V2 == "Refined"])), data = violin.plot.dat,
#          method = 'swarm',
#          pch = 16,
#          xlab = '', ylab = 'Follow-up time (months)')#,
# #         labels = c('Censored', 'Metastasis'))

# Trying to create a scatterplot of the differences in predicted (MSPminer vs refined SG) and the true abundance

plot(colSums(init_tax_abund), colSums(true_tax),xlim = c(0, 450), ylim = c(0,250))
abline(a=0, b=1) # adding a x=y line
points(colSums(final_tax_abund), colSums(true_tax),col="green")
#points(colSums(true_tax), colSums(true_tax),col="red")



# Creating a bubble plot to visualize the relative abundances 


# converting the matrix to a df

init.read.df <- data.frame(x=rep((dimnames(init_tax_abund)[[1]]),each=ncol(init_tax_abund)),y=rep((dimnames(init_tax_abund)[[2]]),times=nrow(init_tax_abund)),z=as.vector(init_tax_abund))
final.read.df <- data.frame(x=rep((dimnames(final_tax_abund)[[1]]),each=ncol(final_tax_abund)),y=rep((dimnames(final_tax_abund)[[2]]),times=nrow(final_tax_abund)),z=as.vector(final_tax_abund))
true_tax.df <- data.frame(x=rep((dimnames(true_tax)[[1]]),each=ncol(true_tax)),y=rep((dimnames(true_tax)[[2]]),times=nrow(true_tax)),z=as.vector(true_tax))


par(mfrow = c(1,2))#,xpd=TRUE)
p1 <- ggplot(data=init.read.df,aes(x=x,y=y,fill=z))+
  geom_tile()+theme(axis.text.y=element_blank(), legend.position = "none") #+ scale_fill_gradient(low="white", high="blue")
p2 <- ggplot(data=final.read.df,aes(x=x,y=y,fill=z))+
  geom_tile()+theme(axis.text.y=element_blank(),legend.position = "none") #+ scale_fill_gradient(low="white", high="blue")
p3 <- ggplot(data=true_tax.df,aes(x=x,y=y,fill=z))+
  geom_tile()+theme(axis.text.y=element_blank(),legend.position = "none") #+ scale_fill_gradient(low="white", high="blue")

grid.arrange(p1, p2, p3, ncol=3)

ggplot(data = init.read.matrix.df %>% 
         summarize(z = n(), AvgLevel = mean(as.integer(z))),
       aes(x = init.read.matrix.df$x, y = init.read.matrix.df$y, size = init.read.matrix.df$z, col = AvgLevel)) +
  geom_point() + theme_bw(base_size = 18) # +
  # scale_colour_gradientn(
  #   colours  = rev(topo.colors(2)),
  #   na.value = "transparent",
  #   breaks   = as.integer(MusicianInterestsSmall$Level) %>% 
  #     unique %>% sort,
  #   labels   = levels(MusicianInterestsSmall$Level),
  #   limits   = c(as.integer(MusicianInterestsSmall$Level) %>% min,
  #                as.integer(MusicianInterestsSmall$Level) %>% max)) +
  # scale_size_continuous(range = c(3, 11)) 
# 
# library(metacoder)
# 
# 
# #load ggplot2
# library(ggplot2)
# df <- data.frame(MSPminer=as.vector(init_true), refined=as.vector(final_true), )
# #create new column for average measurement
# df$avg <- rowMeans(df) 
# 
# #create new column for difference in measurements
# df$diff <- df$MSPminer - df$refined
# 
# #view first six rows of data
# head(df)
# 
# #find average difference
# mean_diff <- mean(df$diff)
# 
# mean_diff
# 
# #find lower 95% confidence interval limits
# lower <- mean_diff - 1.96*sd(df$diff)
# 
# lower
# 
# #find upper 95% confidence interval limits
# upper <- mean_diff + 1.96*sd(df$diff)
# 
# upper
# 
# #create Bland-Altman plot
# ggplot(df, aes(x = avg, y = diff)) +
#   geom_point(size=2) +
#   geom_hline(yintercept = mean_diff) +
#   geom_hline(yintercept = lower, color = "red", linetype="dashed") +
#   geom_hline(yintercept = upper, color = "red", linetype="dashed") +
#   ggtitle("Bland-Altman Plot") +
#   ylab("Difference Between Measurements") +
#   xlab("Average Measurement") +  scale_y_continuous(trans='log10')
# 
# plot(true_tax,abs(final_true),col='#0571B0', ylim=c(0,40), log="x") #xlim=c(0,25)
# par(new=TRUE)
# plot(true_tax,abs(init_true),col='#999999', ylim=c(0,40), log="x")


# The figure in the article

df <- data.frame(MSPminer=as.vector(abs(initial_true)), Refined=as.vector(abs(final_true)), Abundance=as.vector(true_tax))
colors <- c("MSPminer" = "black", "Refined" = '#92C5DE')
p1 <- ggplot(data=df,  aes(x=Abundance, y = Refined, color="Refined")) + geom_point() + geom_point(data=df, aes(x=Abundance, y = MSPminer, color="MSPminer")) +
  xlab("True relative abundance (log)") + theme_minimal() +xlim(-3,50) + scale_x_continuous(trans='log10') + 
  ylab("Difference in relative abundance") +# ylim(0,50)+
  ggtitle("Error in relative abundance given the calculated and true abundance of species-level taxa") + #, subtitle = "Each bar represents a sample. Samples are arranged by increasing abundance of 'Ruminococcus sp. SR1/5'.")+
  scale_color_manual(values = colors) + labs(color="Method")
p1



# p2 <- ggplot(data=df,  aes(x=Abundance, y = Refined)) + geom_point() +
#   xlab("True Abundance") + theme_minimal() + xlim(-3,50) +scale_x_continuous(trans='log10') + 
#   ylab("Difference in abundance") + ylim(0,50)+
#   ggtitle("B: Refined Signature Genes") + #, subtitle = "Each bar represents a sample. Samples are arranged by increasing abundance of 'Ruminococcus sp. SR1/5'.")+
#   #  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(size = 11))# +
# #  labs(fill='Species') + scale_y_continuous(labels = function(x) paste0(x*100, "%"),limits = c(0, 1), expand = c(0, 0)) 
# 
# 
# abundance_plots <- ggarrange(p1, p2 + theme(axis.text.y = element_blank(),
#                                             axis.ticks.y = element_blank(),
#                                             axis.title.y = element_blank() ), ncol=2, nrow=1, widths=c(1,0.91), common.legend=TRUE, legend = "bottom")
# 
# annotate_figure(abundance_plots, top = text_grob("Error in relative abundance given the true abundance of species-level taxa", 
#                                                  color = "black", face = "bold", size = 15))





# Testing the abundance using all genes from MSPminer (not onlky the SG/initial 100 genes)

mspminer.read.matrix <- matrix(NA, nrow=dim(Clusterlist[[1]])[2], ncol=length(Clusterlist))
sample.names <- colnames(Clusterlist[[1]])
rownames(mspminer.read.matrix) <- sample.names # setting the sample names as rownames
colnames(mspminer.read.matrix) <- names(Clusterlist) # setting the cluster ID's as colnames
mspminer.Clusterlist <- Clusterlist



for (id in names(Clusterlist)){
  #removing the samples, with less than 3 genes with mapped reads
  #init.gene.names <- outputs_parallel[[id]]$genes$init
  init.colsum <- colSums(Clusterlist[[id]][, ])
  init.mapped <- names(init.colsum[init.colsum >= n.mapped.minimum])
  init.not_mapped <- sample.names[!sample.names %in% init.mapped]
  mspminer.Clusterlist[[id]][,init.not_mapped] <- 0 # setting the counts to 0 if less than n.mapped.minimum genes have reads that map
  
  # The readcounts are divided by the gene length
  init.reads <- mspminer.Clusterlist[[id]][, ] / GeneLengths[rownames(mspminer.Clusterlist[[id]])]
  
  # summing the read counts for the id/cluster/MGS
  mspminer.read.matrix[, id] <- colSums(init.reads)
}

mspminer.abundance <- mspminer.read.matrix/rowSums(mspminer.read.matrix)
mspminer_tax_abund <- matrix(0, nrow = dim(true_tax)[1], ncol = dim(true_tax)[2], dimnames=list(row.names(true_tax), colnames(true_tax)))

for (j in 1:length(colnames(true_tax))){
  col = colnames(true_tax)[j]
  col = stringr::str_replace_all(col, "\\s", " ")  #changing the encoding
  for (i in 1:length(tax_df$Taxonomy)){
    if (col == tax_df$Taxonomy[i]){
      mspminer_tax_abund[,j] <- mspminer_tax_abund[j] + mspminer.read.matrix[,tax_df$MSP[i]]
    }
  }
}


# Normalizing
mspminer_tax_abund_norm <- ((mspminer_tax_abund)/rowSums(mspminer_tax_abund)*100)

# Subtracting the true abundance from the predicted 
mspminer_true <- mspminer_tax_abund_norm-true_tax

# Testing whether the predicted relative abundance distributions are identical
wilcox.test(mspminer_true, final_true, paired = TRUE, conf.int=TRUE)
wilcox_effsize(data.frame("init_true"=as.vector(init_true), "final_true"=as.vector(final_true)), init_true ~ final_true) # takes a long time to run


df <- data.frame(MSPminer=as.vector(abs(init_true)), Refined=as.vector(abs(final_true)), Abundance=as.vector(true_tax)) # , All_MSPminer=as.vector(abs(mspminer_true))
colors <- c("MSPminer" = "black", "Refined" = '#92C5DE') #, All_MSPminer="red")
p1 <- ggplot(data=df,  aes(x=Abundance, y = Refined, color="Refined")) + geom_point(size=1) + geom_point(data=df, aes(x=Abundance, y = MSPminer, color="MSPminer"),size=1) +# geom_point(data=df, aes(x=Abundance, y = All_MSPminer, color="All_MSPminer")) +
  xlab("True relative abundance (log)") + theme_minimal() +xlim(-3,50) + scale_x_continuous(trans='log10') + 
  ylab("Difference in relative abundance") + # ylim(0,30)+
  ggtitle("Error in relative abundance given the calculated and true abundance of species-level taxa") +
  scale_color_manual(values = colors) + labs(color="Method") 
p1
