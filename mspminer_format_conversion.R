# Takes in the read counts of the genes, MSPminer's bin-assignment of the genes and the depth file computed by the jgi_summarize_bam_contig_depth script provided by MetaBAT
# All found at https://zenodo.org/record/4306051#.YgpJLC8w2Al and provided by: 
# Borderes M, et al. A comprehensive evaluation of binning methods to recover human gut microbial species from a non-redundant reference gene catalog. NAR Genom Bioinform. 2021 Mar 1;3(1):lqab009.

# Loading the genes and their bins predicted by MSPminer
all_MSPs <- read.table("MSPminer.binning", header = F, sep = "\t", skip=4)
colnames(all_MSPs) <- c("SEQUENCEID", "MSPID")

# Load in the read counts
genemat <- read.table("SGC_base_norm.tsv", header = T, sep = "\t", row.names=1) # base normalized read counts
rows <- length(genemat[,1])


#bincounts <- as.data.frame(table(genes_bins[,3])) #counting the genes of all bins/samples!


occurences <- table(unlist(all_MSPs$MSPID)) # getting the number of genes in each MSP

# For every MSP, assign the genes and it's respective counts and order according to coabundance
for (i in 1:length(occurences)){
  gene_names <- all_MSPs[all_MSPs$MSPID == names(occurences[i]),1] # getting the gene names of the MSP
  
  # Sort the gene count matrix of the MSP according to correlation
  AbuMatrix <- as.matrix(genemat[gene_names,])
  Median_profile <- colMedians(AbuMatrix)
  pcc2MedianProfile <-cor(t(AbuMatrix) , Median_profile )
  gene_order <- row.names(pcc2MedianProfile)[order(pcc2MedianProfile, decreasing = TRUE)]
  ordered_AbuMatrix <- AbuMatrix[gene_order,]
  assign(names(occurences[i]), ordered_AbuMatrix)
}

#combining the MSPs into a large list and save the object
MSPlist <- mget(names(occurences))
saveRDS(MSPlist, file="MSP_coabundance_SGC.RDS")


# Saving the gene lengths as a vector
gene_lengths <- read.table("SGC_depth_jgi.tsv", header=T, sep="\t")
genes <- split(gene_lengths$contigLen, gene_lengths$contigName)
genes <- unlist(genes)

saveRDS(genes, file="SGC_genelengths.RDS")
