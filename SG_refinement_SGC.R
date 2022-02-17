# Takes in the MSPobject, the genelengths and the functions as input.
# 
# Identifies an optimal Signature gene (SG) set and outputs both a new MSPobject as well as a csv with the MSE of the initial and refined SG set. 


# initialization and parameter-tuning
library(ggplot2) 
library(vcd) #goodfit
library(plyr)
library(VGAM) # for gpois
library(matrixStats) #rowMedians
library(parallel) #for paralleling; mclapply

set.seed(1337)
window_size <- 699
turkey_factor <- 1.2

work_dir = "/Users/trinezachariasen/Documents/PhD-Metagenomics/MSPminer_benchmark/Github_versions/"
setwd(work_dir)

GeneLengths <- readRDS('SGC_genelengths.RDS')
MSPlist <- readRDS("MSP_coabundance_SGC.RDS") 
source("Functions.R")

#Number of signature genes
n.genes <- 100

#minimum number of mapped genes required
n.mapped.minimum <- 3

ids <- names(MSPlist) #the ids of the MGS

# initializing objects to keep track of stats throughout the process
mse.before <- rep(0, length(ids)) # MSE of the initial signature gene set according to the expected distribution
names(mse.before)  <- ids
mse.after <- (mse.before)  # MSE of the refined signature gene set
mse.final <- (mse.before)  # MSE of the refined signature gene set after removal of outliers


# looping over the ids of the MGS
#registers time to run each
mapped_samples <- c()
n_samples <- c()

# identifying the names of genes in all MGS the data set
MSP_gene_names <- c()
for(MSP in names(MSPlist)){
  MSP_gene_names <- c(MSP_gene_names, rownames(MSPlist[[MSP]]))
}

present_genes <- GeneLengths[names(GeneLengths) %in% MSP_gene_names]

# counting the number of genes of each MSP 
gene_count <- c()
for (MSP in names(MSPlist)){
     gene_count <- c(gene_count, length(MSPlist[[MSP]])/ncol(MSPlist[[1]]))
     if ((length(MSPlist[[MSP]])/ncol(MSPlist[[1]]))<100){print(MSP)}
}


# Prescreening the genes, sorting them according to be uniformly identified
if(file.exists(paste(work_dir, "trimmed_MSP_SGC.RData", sep = ""))){
  load(paste(work_dir, "trimmed_MSP_SGC.RData", sep = ""))
   MSPlist <- MSPlist_trimmed
   rm(MSPlist_trimmed)
}else{
  #preselecting genes
  print('pre-screeing genes')
  for(id in ids){

    genes_r <- MSPlist[[id]]
    gene_count <- length(MSPlist[[id]])/ncol(MSPlist[[1]])
    
    if (gene_count>100){ # we want at LEAST 100 unique genes
      final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))
  
      #finds the frequency of each gene. We want a ´p_sig_gene~1/n_signature_gene
      #Because if we just have one that's messy and  more likely than the others, then you effed. having a high one is actually worse, since 95 or more is fine
      frequency_of_genes <- rowSums(final.reads>0, na.rm = T)
  
      #So! Basically...
      #if I don't exclude certain genes, we end up with scenarios in which you get a 3-4 genes that are 5x as likely to be identified compared to the ones that are hard to find
      #That means, that a LOT of instances where if you sequence, say, 10 times, anything below nine is significant given a totally random draw. But!
      #If you have 3-4 genes that are very likely to pop up, you end up getting a ton of times where the same high-frequency gene gets found TWICE!
      #So, we want to choose a set of genes we can purify on that 1) are usually found in high frequency and 2) are relatively consistent with each other.
  
      #in a test, I tried to select all genes between median*1.2 and median/1.2 and only run on those. Bad sutff. a lot of high-read samples never got all genes, good stuff: they followed a n_expected_signature_genes =80 extremely well
      #But, of course, I did choose a lot of genes that are sort of middling in terms of frequency, i.e. relatively hard to detect, indicating that they probably aren't found in all members of that MGS...
      #So, let's try to only select high-frequency genes, none of which are outliers by turkey's method, but a bit more sensitive.
      vars <- c()
      good_set_found <- F
      if(length(frequency_of_genes)>window_size+1){
        frequency_of_genes <- sort(frequency_of_genes,decreasing = T)
        for(i in seq(1,length(frequency_of_genes)-window_size)){
          in_window <- frequency_of_genes[seq(i,i+window_size)]
          iqr <- IQR(in_window,na.rm = T)
          upper_bound <- median(in_window,na.rm = T)+iqr*turkey_factor
          lower_bound <- median(in_window,na.rm = T)-iqr*turkey_factor
          
          vars <- c(vars,var(in_window),na.rm=T)
          #fewer than 1% of samples are outliers
          if(sum(in_window>upper_bound | in_window<lower_bound)<=window_size/100){
            good_set_found <- T
            break
          }
        }
  
        if(good_set_found==F){
          names(vars) <- seq(1,length(frequency_of_genes)-window_size)
          best_vars <- na.omit(names(vars)[vars==min(vars)])
          if(length(best_vars)>0){
            i <- best_vars[1]
          }
        }
      }
      #if you can't find a good region, go on then
      if(good_set_found==T){
        acceptable_genes <- names(frequency_of_genes)[seq(i,window_size+i)]
  
        original_genes <- rownames(MSPlist[[id]])[1:100]
  
        #adds the original ones. they'r good for SOMETHING and reorders them to be the original order just to be safe
        acceptable_genes <- unique(c(original_genes,acceptable_genes))
        acceptable_genes <- rownames(MSPlist[[id]])[rownames(MSPlist[[id]])%in%acceptable_genes]
        #trims the set of possible genes down to acceptable ones
        MSPlist[[id]] <- MSPlist[[id]][acceptable_genes,]
  
      }
    } else{
      MSPlist[[id]] <- NULL # if there are <100 unique gene in the entity
    }
    

  }
  MSPlist_trimmed <- MSPlist
  save(MSPlist_trimmed, file = paste(work_dir, "trimmed_MSP_SGC.RData", sep = ""))
}


ids <- names(MSPlist)

numCores <- 2 # Number of threads used for the gene optimization

# Running the actual Signature Gene refinement of the ids using the function run_one_id
outputs_parallel <- mclapply(ids, run_one_id, mc.cores = numCores)
names(outputs_parallel) <- ids

# Storing the output as a dataframe
mgs_results = ids
names(outputs_parallel) <- mgs_results

MSE_df = data.frame(matrix(unlist(outputs_parallel)[grep("*.MSE.", names(unlist(outputs_parallel)))], nrow=length(mgs_results), byrow=T))
row.names(MSE_df) <- names(outputs_parallel)
colnames(MSE_df) <- c("Initial gene set", "Intermediate", "Final gene set")
print(MSE_df)
write.csv(MSE_df, 'mse.csv') # Writing the results to a csv
 
# storing the biological entities along with its genes
saveRDS(outputs_parallel, file="MSP_SGC_refined_coabundance.RDS") 
