plot_general_illustration=T

plot_all_poissons=T
#in the article I hope to have a few cases that are illustrative of how the improvements are made
MGSs_of_interest <- list()

MGSs_of_interest[['big_improvement_big_change']]<- list()
MGSs_of_interest[['big_improvement_little_change']]<- list()
MGSs_of_interest[['little_improvement_big_change']]<- list()
MGSs_of_interest[['just_ignore_the_bad_samples']]<- list()
MGSs_of_interest[['big_change_MSE']]<- list()

n_signature_genes_expected <- '95'
minimum_sampels <- 3

#presets
set1 <- 'init'
set2 <- 'best'

#sankey stuff
library(networkD3)
library(htmlwidgets)
library(htmltools)
library(gridGraphics)
library(cowplot) #ggdraw
library(ggplot2)
library(ggrepel)

#loading CM colors
cmcol1 = c('#0571B0', '#92C5DE', '#999999', '#E0E0E0', '#CA0020',
           '#FF8784', '#FFD900', '#FFE967', '#349F2C', '#BAE4B3',
           '#6B3D9A', '#CAB2D6')
#sum(coupon_collector_simulations[[100]][[85]][names(u_table)<50])/sum(coupon_collector_simulations[[100]][[85]])
signature_genes_test_simulate <- function(n_signature_genes,n_reads,perm=10^6){
  urn = seq(1:n_signature_genes)
  u = replicate(perm, length(unique(sample(urn,n_reads,repl=T))))
  return(u)
}

#capitalizes string
proper=function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

coupon_sim_pval_just_analytical<- function(n_reads,n_signature_genes_expected,observed_signature_genes,mode){
  #gets the p-value of the amount of observed signature genes given a certain number of reads assigned to signature genes or fewer/greater
  #n_reads is the number of reads assigned to signature genes, an int
  #n_signature_genes_expected the amount of signature genes of which we are willing to accept an MGS as being present (ideally 100/100 would be expected, but if it has 99, will we say it DEFINETLY doesn't exist??)
  #observed_signature_genes number of signature genes we observe
  #wether we want to get the p-value for observed_signature_genes or fewer or greater. Options "fewer" and "greater". Default is "fewer"
  #coupon_collector_simulations - the coupon_collector_simulations element. Set by default
  
  #if there are no reads, there is no p-value. Duh
  if(n_reads==0){
    return(NaN)
  }
  
  #first it tries to get an analytical answer. That would always be nice
  returnme <- NaN
  if(n_reads<=dim(analytical_coupon_pvals_array)[[2]]){
    if(mode=='exact'){
      returnme <- analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(observed_signature_genes)]
    }else if(mode=='fewer'){
      returnme <- sum(analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(seq(1,observed_signature_genes,1))])
    }else if(mode=='greater'){
      returnme <- sum(analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(seq(observed_signature_genes,dim(analytical_coupon_pvals_array)[[3]],1))])
      
    }
  }
  return(returnme)
}

coupon_sim_pval_analytical_if_possible <- function(n_reads,n_signature_genes_expected,observed_signature_genes,mode){
  #gets the p-value of the amount of observed signature genes given a certain number of reads assigned to signature genes or fewer/greater
  #n_reads is the number of reads assigned to signature genes, an int
  #n_signature_genes_expected the amount of signature genes of which we are willing to accept an MGS as being present (ideally 100/100 would be expected, but if it has 99, will we say it DEFINETLY doesn't exist??)
  #observed_signature_genes number of signature genes we observe
  #wether we want to get the p-value for observed_signature_genes or fewer or greater. Options "fewer" and "greater". Default is "fewer"
  #coupon_collector_simulations - the coupon_collector_simulations element. Set by default
  
  #if there are no reads, there is no p-value. Duh
  if(n_reads<=0){
    return(NaN)
  }
  
  max_k_straps <- max(as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]])))
  #if you exceed the number of boostrapped k-values, return 0 for n==d and 1 for n!=d
  if(n_reads>=max_k_straps){
    
    if(as.numeric(n_signature_genes_expected)<=observed_signature_genes){
      return(1)
    }else{
      return(0)
    }
  }
  
  
  
  #first it tries to get an analytical answer. That would always be nice
  returnme <- NaN
  if(n_reads<=dim(analytical_coupon_pvals_array)[[2]]){
    if(mode=='exact'){
      returnme <- analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(observed_signature_genes)]
    }else if(mode=='fewer'){
      returnme <- sum(analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(seq(1,observed_signature_genes,1))])
    }else if(mode=='greater'){
      returnme <- sum(analytical_coupon_pvals_array[n_signature_genes_expected,as.character(n_reads),as.character(seq(observed_signature_genes,dim(analytical_coupon_pvals_array)[[3]],1))])
      
    }
  }
  
  if(is.na(returnme)){
    #if you don't have the simulation you need and it's a reasonable request (i.e. not 0 and not above the maximum) then run that simulation real quick
    if(n_reads < max(as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]]))) &
       n_reads!=0 & 
       !n_reads%in%as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]]))){
      message('Oops! Found an missing read number in Coupon collector simulations. Read number: ',as.character(n_reads),'... Simulating now at a lower resolution...')
      temp <- table(signature_genes_test_simulate(n_signature_genes = n_signature_genes_expected,n_reads = n_reads,perm = 10^4))
    }else{
      temp <- coupon_collector_simulations[[as.character(n_signature_genes_expected)]][[as.character(n_reads)]]
    }
    
    
    temp_names <- as.numeric(names(temp))
    
    if(mode=='fewer'){
      return(sum(temp[as.character(temp_names[temp_names<=observed_signature_genes])])/sum(temp))
    }else if(mode=='greater'){
      return(sum(temp[as.character(temp_names[temp_names>=observed_signature_genes])])/sum(temp))
    }else if(mode=='exact'){
      return(sum(temp[as.character(temp_names[temp_names==observed_signature_genes])])/sum(temp))
    }else{
      return(NA)
    }
  }else{
    return(returnme)
  }
  
}

coupon_sim_pval_just_bootstrap <- function(n_reads,n_signature_genes_expected,observed_signature_genes,mode){
  #gets the p-value of the amount of observed signature genes given a certain number of reads assigned to signature genes or fewer/greater
  #n_reads is the number of reads assigned to signature genes, an int
  #n_signature_genes_expected the amount of signature genes of which we are willing to accept an MGS as being present (ideally 100/100 would be expected, but if it has 99, will we say it DEFINETLY doesn't exist??)
  #observed_signature_genes number of signature genes we observe
  #wether we want to get the p-value for observed_signature_genes or fewer or greater. Options "fewer" and "greater". Default is "fewer"
  #coupon_collector_simulations - the coupon_collector_simulations element. Set by default
  
  #if there are no reads, there is no p-value. Duh
  if(n_reads<=0){
    return(NaN)
  }
  
  max_k_straps <- max(as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]])))
  #if you exceed the number of boostrapped k-values, return 0 for n==d and 1 for n!=d
  if(n_reads>max_k_straps){
    if(as.numeric(n_signature_genes_expected)<=observed_signature_genes){
      return(1)
    }else{
      return(0)
    }
  }
  
  
  #if you don't have the simulation you need and it's a reasonable request (i.e. not 0 and not above the maximum) then run that simulation real quick
  if(n_reads < max(as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]]))) &
     n_reads!=0 & 
     !n_reads%in%as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]]))){
    message('Oops! Found an missing read number in Coupon collector simulations. Read number: ',as.character(n_reads),'... Simulating now at a lower resolution...')
    temp <- table(signature_genes_test_simulate(n_signature_genes = n_signature_genes_expected,n_reads = n_reads,perm = 10^4))
  }else{
    temp <- coupon_collector_simulations[[as.character(n_signature_genes_expected)]][[as.character(n_reads)]]
    
  }
  
  temp_names <- as.numeric(names(temp))
  
  if(mode=='fewer'){
    return(sum(temp[as.character(temp_names[temp_names<=observed_signature_genes])])/sum(temp))
  }else if(mode=='greater'){
    return(sum(temp[as.character(temp_names[temp_names>=observed_signature_genes])])/sum(temp))
  }else if(mode=='exact'){
    return(sum(temp[as.character(temp_names[temp_names==observed_signature_genes])])/sum(temp))
  }else{
    return(NA)
  }
}


pois.plot <- function(Genes, Reads, mse,id){ #},HGMGS){
  # Plotting the reads and mapped genes for each sample
  #
  # Args:
  #   Genes: A named vector containing the number of genes mapped in the samples
  #   Reads: Vector containing the number of reads mapping to the Genes. Genes and Reads must 
  #           have the same length, greater than one, with no missing values
  #   mse:  The MSE of the MGS when compared to the expected distribution
  #   id is the ID of the species in the HGMGS object
  #   HGMGS is the MGS object in use
  #
  # Returns:
  #   A plot of all samples with their corresponding genes and mapped reads (log)
  
  
  df_plot <- do.call(rbind, Map(data.frame, Reads = Reads, Genes = Genes, name = names(Genes)))
  p = ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
             data = df_plot, aes(text = name, y = Genes, x = Reads)) + 
    geom_point(size = 0.7, aes(color = "Samples")) +
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes, color="Predicted"), size=1, alpha=0.5) + #t he expected distribution
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected') +
    
    ggtitle(paste(taxonomy$MainTax[taxonomy$Cluster==id], id))
#    ggtitle(paste(HGMGS$Tax[id, "MainTax"], id))
  
  plotly.data <- (p + scale_x_log10() + scale_colour_manual(values = c("Samples" = cmcol1[2], "Outlier" = cmcol1[5], "Predicted"="black")))
  plotly.data <- plotly.data + annotate("text", x=2, y=80, label = paste("MSE:", round(mse, 2), parse=T)) +
  plotly.data <- plotly.data  + theme_minimal() + theme(legend.position="bottom")
  return(plotly.data)  # for creating PDF (not interactive)
}

pois.plot.with.probabilities <- function(Genes, Reads, mse,id,HGMGS,mode,n_signature_genes_expected){
  # Plotting the reads and mapped genes for each sample
  #
  # Args:
  #   Genes: A named vector containing the number of genes mapped in the samples
  #   Reads: Vector containing the number of reads mapping to the Genes. Genes and Reads must 
  #           have the same length, greater than one, with no missing values
  #   mse:  The MSE of the MGS when compared to the expected distribution
  #   id is the ID of the species in the HGMGS object
  #   HGMGS is the MGS object in use
  #   mode is the p-value for "fewer" or "greater"
  #   n_signature_genes_expected number of ints we are willing to accept as "the same"
  # Returns:
  #   A plot of all samples with their corresponding genes and mapped reads (log)
  
  #if you have more reads than this, then the probability of getting anything other than only n_expected is zero, so why bother
  
  
  #gets the p-values
  #gets the p-values
  p_values <- c()
  for(sample in names(Reads)){
    p_values <- c(p_values,
                  CCP_pval(k = Reads[sample],
                          n = n_signature_genes_expected,
                          d = Genes[sample],mode = mode))
  }
  names(p_values) <- names(Reads)
  
  df_plot <- do.call(rbind, Map(data.frame, Reads = Reads, Genes = Genes,pval=p_values, name = names(Genes)))
  
  p = ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
             data = df_plot, aes(text = name, y = Genes, x = Reads)) + 
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes), size=1, alpha=0.2,colour="#b8dde3") + #t he expected distribution
    geom_point(size = 0.7, aes(color = pval)) +
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected') +
    #ggtitle(paste(HGMGS$Tax[id, "MainTax"], id))
    ggtitle(paste(taxonomy$MainTax[taxonomy$Cluster==id],id))
  
  p <- (p + scale_x_log10())
  p <- p + annotate("text", x=2, y=80, label = paste("MSE:", mse), parse=T)
  p <- p  + theme_minimal() + theme(legend.position="bottom")+scale_colour_gradientn("Chance of complete\n signature gene set",colours = wesanderson::wes_palette("Darjeeling1", 1000, type = "continuous"))
  
  
  return(p)  # for creating PDF (not interactive)
  #print(ggplotly(plotly.data, tooltip = c("x", "y", "text")))  #creating interactive plot with plotly
}

#returns plotdata instead of plot for MAXIMUM COSTUMIZATION
pois.plotdata.with.probabilities <- function(Genes, Reads, mse,id,HGMGS,mode,n_signature_genes_expected){
  # Plotting the reads and mapped genes for each sample
  #
  # Args:
  #   Genes: A named vector containing the number of genes mapped in the samples
  #   Reads: Vector containing the number of reads mapping to the Genes. Genes and Reads must 
  #           have the same length, greater than one, with no missing values
  #   mse:  The MSE of the MGS when compared to the expected distribution
  #   id is the ID of the species in the HGMGS object
  #   HGMGS is the MGS object in use
  #   mode is the p-value for "fewer" or "greater"
  #   n_signature_genes_expected number of ints we are willing to accept as "the same"
  # Returns:
  #   A plot of all samples with their corresponding genes and mapped reads (log)
  
  #if you have more reads than this, then the probability of getting anything other than only n_expected is zero, so why bother
  maximum_reads_allowed <- max(as.numeric(names(coupon_collector_simulations[[as.character(n_signature_genes_expected)]])))
  
  
  #gets the p-values
  p_values <- c()
  for(sample in names(Reads)){
    p_values <- c(p_values,
                  CCP_pval(k = Reads[sample],
                        n = n_signature_genes_expected,
                        d = Genes[sample],mode = mode))
  }
  names(p_values) <- names(Reads)
  
  df_plot <- do.call(rbind, Map(data.frame, Reads = Reads, Genes = Genes,pval=p_values, name = names(Genes)))
  
  return(df_plot)  # for creating PDF (not interactive)
  #print(ggplotly(plotly.data, tooltip = c("x", "y", "text")))  #creating interactive plot with plotly
}

gene.selection <-  function(n.replace, gene.prob, gene.pos) {
  #  Identifying the best suited replacement signature genes
  # Args:
  #   n.replace: number of genes to replace
  #   gene.prob: Vector with the probability and name of the genes
  #   gene.pos: Vector with the position of the gene
  # Returns:
  #   replace: Vector with genes with a high Pearson residual in one or more samples. 
  
  if (n.replace > length(gene.prob)){
    n.replace <- length(gene.prob)
  }
  
  #replace.genes <- sort(abs(gene.prob) * 0.95 - gene.pos * 0.05, decreasing = TRUE)[0:n.replace]
  replace.genes <- sort(abs(gene.prob) * 0.975 + gene.pos * 0.025, decreasing = FALSE)[0:n.replace]
  
  return("replace" = replace.genes)
}

rank.mat <- function(id,gene.names, step, phase, threshold, t.start, t.end, n.replace = FALSE,mapped) {
  #  Identify suitable signature genes by fitting NB model and ranking the genes
  #
  # Args:
  #   gene.names: A list of genes used to train the model
  #   step: Indicating whether the suitable genes are found based on "median" or "quantile"
  #   phase: The phase of the modelling. Currently there are two phases: "initial" and "rotation"
  #ØP: Soooo... What's the difference?
  #   threshold: Mean value of rank for keeping genes
  #   t.start: The start of the genes in the read count matrix, default: 0 for initial, 100 for rotation
  #   t.end: The end of the genes in the read count matrix, default: 100 for initial, 500 for rotation
  #   n.replacee: The number of genes to replace, int
  #   mapped are the samples we are confident the MGS is present in
  #
  # Returns:
  #   rank.mat: Matrix containing the ranks of the gene.names in the samples
  #   good.genes: List of the gene names of the genes that we keep in the signature gene set, for rotation it is the list of genes of length n.replace
  #   countReads: double with sample name and corresponding total read count
  #   countGenes: double with sample name and corresponding detected number of signature genes
  #   mse: The MSE of the fit to the expected distribution (s=(1-((G-1)/G)^N)*G)
  
  
  #ØP: what is this mapped business? It's not an input?
  rank.matrix <- matrix(NA, nrow=(t.end-t.start), ncol=length(mapped))
  
  reads <- round(utuseq_MGSgeneCountsML[[id]][gene.names, names(mapped)] / 
                   (present_genes[rownames(utuseq_MGSgeneCountsML[[id]][gene.names, ])] * 10^-3))
  
  # running goodfit on all samples, saving the model in the gf matrix
  gf.nbin <- apply(matrix(1:length(mapped), nrow = length(mapped), ncol = 1), 1, function(x) (goodfit(reads[, x], type = "nbinom", method = "ML", par = list(size = mean(reads[,x])))))
  
  
  # filling in the rank.matrix with the rank of all genes within each sample
  if (phase == "initial") {
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[reads[,c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
    
  } else if (phase == "rotation") {
    new.reads <- matrix(NA, nrow = (t.end-t.start), ncol = length(mapped))
    new.reads <- round(utuseq_MGSgeneCountsML[[id]][(t.start + 1):t.end, 1:length(names(mapped))] / 
                         (present_genes[rownames(utuseq_MGSgeneCountsML[[id]][(t.start + 1):t.end, ])] * 10^-3))
    
    # if a new gene have more reads than found in the good genes, then set the readcount to the max observed in the good genes
    for (i in 1:length(mapped)) new.reads[, i][new.reads[, i] > max(gf.nbin[[i]]$count)] = max(gf.nbin[[i]]$count)
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[new.reads[, c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
  }
  
  rownames(rank.matrix) <- names(utuseq_MGSgeneCountsML[[id]][(t.start + 1):t.end, 1]) # genes as rownames
  colnames(rank.matrix) <- names(mapped) # samples as colnames
  
  
  if (phase == "initial") {
    
    # identifying the genes to replace
    if (step == "mean"){
      gene.performance <- rowMeans(rank.matrix) # identifying the means of the signature genes
    } else if(step=="quantile"){
      gene.performance <- rowQuantiles(rank.matrix, probs = seq(from = 0, to = 1, by = 0.05))[,20] #identifying the 95-percentile of the signature genes
    }
    
    low.prob.genes <- gene.names[gene.performance > threshold] # if the rank of the genes are above threshold
    n.replace <- length(low.prob.genes) # number of genes to replace
    good.genes <- gene.names[!(gene.names %in% low.prob.genes)] # the genes with rank below the threshold
    
    
    countReads <- colSums(round((utuseq_MGSgeneCountsML[[id]][1:n.genes, ] /
                                   (present_genes[rownames(utuseq_MGSgeneCountsML[[id]][1:n.genes, ])] * 10^-3))));
    
    countGenes <- colSums(utuseq_MGSgeneCountsML[[id]][1:n.genes, ] > 0)
    
  } else if (phase == "rotation"){
    if (step == "mean"){
      gene.performance <- rowMeans(rank.matrix)
    } else if(step=="quantile"){
      gene.performance <- rowQuantiles(rank.matrix, probs = seq(from = 0, to = 1, by = 0.05))[,20] #95percentile
    }
    names(gene.performance) = names(utuseq_MGSgeneCountsML[[id]][(t.start + 1):t.end, 1]) 
    high.rank <- gene.performance 
    high.pos <- seq((t.start + 1), t.end)
    
    # ensuring that the new genes are not already part of the signature gene set
    high.rank <- high.rank[! names(gene.performance) %in% gene.names]
    high.pos <- high.pos[! names(gene.performance) %in% gene.names]
    
    # identifying the best replacement genes
    good.genes <- gene.selection(n.replace, high.rank, high.pos)
    
    replace.reads <- utuseq_MGSgeneCountsML[[id]][names(good.genes), , drop = FALSE]
    genes_r <- rbind(utuseq_MGSgeneCountsML[[id]][gene.names, ], replace.reads)
    
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))
    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
  }
  
  # calculating the predicted genecounts
  pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
  
  # calculating the MSE
  mse <- round(mean((countGenes - pred)^2), 2)
  
  return(list("good.genes" = good.genes, "countReads" = countReads, "countGenes" = countGenes, "mse" = mse))
}

plot_relevant_pois <- function(MGS){
  outlist <- list()
  
  entries_of_interest <- names(trine_MGS_object[['i']])[grepl(MGS,names(trine_MGS_object[['i']]))]
  entries_of_interest <- entries_of_interest[!grepl('random',entries_of_interest)]
  for(entry in entries_of_interest){
    genes_of_interest <- names(present_genes[trine_MGS_object[["i"]][[entry]][1:n.genes]])
    colsum <- colSums(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest,])
    mapped <- colsum[colsum >= 3]
    
    
    reads <- round(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, names(mapped)] / 
                     (present_genes[rownames(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ])] * 10^-3))
    
    #gets the total number of reads for each sample
    countReads <- colSums(round((utuseq_MGSgeneCountsML[[MGS]][genes_of_interest,] /
                                   (present_genes[rownames(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ])] * 10^-3))));
    #getting the number of identified signature genes for each sample
    countGenes <- colSums(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ] > 0)
    
    # calculating the predicted genecounts and then MSE
    pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
    
    # calculating the MSE
    mse <- round(mean((countGenes - pred)^2), 2)
    
    #pois.plot(init.fit$countGenes, init.fit$countReads, init.fit$mse,MGS,HGMGS) # The initial signature gene set
    outlist[[entry]] <- pois.plot(Genes = countGenes,
                                  Reads = countReads,
                                  mse = mse,
                                  id = MGS,
                                  HGMGS = HGMGS)+
      labs(subtitle=gsub('mean','Mean-rank refined',gsub('best','best performing genes removed',gsub('init','Initial fit',rev(strsplit(entry,'_')[[1]])[1]))))
  }
  
  return(outlist)
}

plot_relevant_pois_probability <- function(MGS){
  outlist <- list()
  
  entries_of_interest <- names(trine_MGS_object[['i']])[grepl(MGS,names(trine_MGS_object[['i']]))]
  entries_of_interest <- entries_of_interest[!grepl('random',entries_of_interest)]
  for(entry in entries_of_interest){
    genes_of_interest <- names(present_genes[trine_MGS_object[["i"]][[entry]][1:n.genes]])
    colsum <- colSums(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest,])
    mapped <- colsum[colsum >= 3]
    
    
    reads <- round(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, names(mapped)] / 
                     (present_genes[rownames(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ])] * 10^-3))
    
    #gets the total number of reads for each sample
    countReads <- colSums(round((utuseq_MGSgeneCountsML[[MGS]][genes_of_interest,] /
                                   (present_genes[rownames(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ])] * 10^-3))));
    #getting the number of identified signature genes for each sample
    countGenes <- colSums(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ] > 0)
    
    # calculating the predicted genecounts and then MSE
    pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
    
    # calculating the MSE
    mse <- round(mean((countGenes - pred)^2), 2)
    
    #pois.plot(init.fit$countGenes, init.fit$countReads, init.fit$mse,MGS,HGMGS) # The initial signature gene set
    outlist[[entry]] <- pois.plot.with.probabilities(Genes = countGenes,
                                                     Reads = countReads,
                                                     mse = mse,
                                                     id = MGS,
                                                     HGMGS = HGMGS,mode = 'fewer',n_signature_genes_expected = n_signature_genes_expected)+
      labs(subtitle=gsub('mean','Mean-rank refined',gsub('best','best performing genes removed',gsub('init','Initial fit',rev(strsplit(entry,'_')[[1]])[1]))))
  }
  
  return(outlist)
}


#loads nessecary stuff
# initializing
fig_dir = "/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/resubmission/"
work_dir = "/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark"
data_dir = "/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/"

#loads relevant data and sets the output directory

n_reads_minimum <- 3
n.genes <- 100


#loads 
setwd(work_dir)

source("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/Github_versions/R_pvalue_getter.R")
source("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/CM/Artikel/R_Anders/Functions_v4.R")
# #load("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/CM/Artikel/R_Anders/coupon_simulations_low_res.Rdata") #coupon simulations
# load("/Users/trinezachariasen/Downloads/coupon_simulations_low_res.Rdata")
# load(paste(data_dir, "trimmed_MSP_SGC.RData", sep="")) # The refined SG for the clusters
# load(paste(work_dir,'/MSP_debug_object','.Rdata',sep=''))
# trine_MGS_object <- MSP_object
# 
# utuseq_MGSgeneCountsML <- readRDS(paste(data_dir, "MSP_SGC_refined_coabundance_normalized.RDS", sep="")) #read count matrices for the Clusters
# GeneLengths <- readRDS(paste(data_dir, "SGC_genelengths.RDS", sep="")) #VAMB gene lengths 
# tax_df <- readRDS("/Users/trizac/Dropbox/PhD-Metagenomics/MSPminer_benchmark/tax_df.RDS")
# #taxonomy <- read.csv("/Users/trizac/Dropbox/CM/VAMB/GTDB-tk_annotated_clusters.csv", header=FALSE) # the taxonomy 
# #colnames(taxonomy) <- c("Cluster", "Taxonomy")

#source("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/Functions_v4_nogenelength.R")
load("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/CM/Artikel/R_Anders/coupon_simulations_low_res.Rdata")
load("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/trimmed_MSP_SGC.RData")
load("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/MSP_debug_object.Rdata")
trine_MGS_object <- MSP_object

utuseq_MGSgeneCountsML <- readRDS("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/PhD-Metagenomics/MSPminer_benchmark/MSP_normalized_coabundance_SGC.RDS") #read count matrices for the Clusters
GeneLengths <- readRDS(paste(data_dir, "SGC_genelengths.RDS", sep="")) #VAMB gene lengths 
tax_df <- readRDS(paste(data_dir, "tax_df.RDS", sep=""))


gene_index <- seq(1,length(GeneLengths))

gene_names <- names(GeneLengths)


if(file.exists(paste0(data_dir,'analytical_pvals_coupon_problem.Rdata'))){
  load(paste0(data_dir,'analytical_pvals_coupon_problem.Rdata'))
}else{
  message('Recreating analytical results from coupon problem solution')
  #calculates the exact chance of getting d different coupons given n total number of unique coupons and k draws
  coupon_d_out_n_k_draws_analytical_stored <- function(d,k,n){
    if(n^k==Inf){
      return(NaN)
    }else if(k<=ncol(S2_matrix)){
      S2_number <- S2_matrix[k,d]
      returnme <- choose(n,d)*S2_number*factorial(d)/n^k
    }else{
      returnme <- NaN
    }
    return(returnme)
  }
  
  second_kind_sterling_numbers <- read.delim(paste0(data_dir,"second_kind_sterling_numbers.tsv"))
  S2_matrix <- matrix(0,nrow = max(second_kind_sterling_numbers[,"n"]),ncol = max(second_kind_sterling_numbers[,"k"]))
  for(row in rownames(second_kind_sterling_numbers)){
    S2_matrix[second_kind_sterling_numbers[row,"n"],second_kind_sterling_numbers[row,"k"]] <- second_kind_sterling_numbers[row,"s2.n.k."]
  }
  
  
  dims <- list(n=seq(1,100),k=seq(1,ncol(S2_matrix)),d=seq(1,nrow(S2_matrix)))
  analytical_coupon_pvals_array <- array(0,dim = lengths(dims),dimnames = dims)
  for(n in dims[["n"]]){
    print(n)
    for(k in dims[["k"]]){
      for(d in dims[["d"]]){
        analytical_coupon_pvals_array[as.character(n),as.character(k),as.character(d)] <- coupon_d_out_n_k_draws_analytical_stored(d,k,n)
      }
    }
  }
  rm(second_kind_sterling_numbers,S2_matrix)
  save(analytical_coupon_pvals_array,file = paste0(data_dir,'analytical_pvals_coupon_problem.Rdata'))
}


message('Loading complete')

#####################################################################################
#########################loading over################################################
#####################################################################################

#setting important variables
#which genes are present
utuseq_names <- c()
for(MGS in names(utuseq_MGSgeneCountsML)){
  utuseq_names <- c(utuseq_names, rownames(utuseq_MGSgeneCountsML[[MGS]]))
}

utuseq_names <- unique(utuseq_names)

present_genes <- GeneLengths[names(GeneLengths) %in% utuseq_names]

names_avail <- names(trine_MGS_object[["i"]])
MGSs <- unique(unlist(lapply(names_avail[!grepl('random',names_avail)],function(i){paste(strsplit(i,'_')[[1]][1],strsplit(i,'_')[[1]][2], sep="_")})))

MGSs_of_interest[['Common_ones']] <- sort(MGSs)[c(1,2,3)]

# inserting NA for the Clusters that do not have a annotation
# Identifying the "Maintax"

# taxonomy$MainTax <- rep("NA", length(taxonomy))
# count <- 0 
# for (MGS in taxonomy$Cluster){
#   if (length(taxonomy$Taxonomy[taxonomy$Cluster==MGS])==0){
#     na_df <- data.frame(MGS, "NA;NA;NA", "NA")
#     colnames(na_df) <- c("Cluster", "Taxonomy", "MainTax")
#     taxonomy <- rbind(taxonomy, na_df)
#     count <- count + 1
#   } else {
#     split_taxonomy <- rev(strsplit(taxonomy$Taxonomy[taxonomy$Cluster==MGS],";")[[1]])
#     main_tax <- c()
#       
#     # as we don't kniw the level of taxonomy we go through each level
#     for (tax in split_taxonomy){
#       
#       #check if the taxonomy is not empty
#       length_tax <- nchar(tax)
#       
#       if (length_tax >= 4){ #if the tax is not empty
#         main_tax <- paste(tax, main_tax)
#         
#         if (is.null(nchar(main_tax))){ # if we have 2 tax in the main tax, then we are done
#           skip }
#         if (nchar(main_tax) > (length_tax+1)){
#             break
#           }
#       }
#     }
#   
#     
#     #main_tax <- paste(second_low_tax, low_tax)
#     taxonomy$MainTax[taxonomy$Cluster == MGS] <- main_tax
#     
#   }
# }




overlap_table <- matrix(NA,ncol=4,nrow=0)
colnames(overlap_table) <- c('MGS','set1','set2','overlap')
#for each MGS get initial an mean performance


#number of samples with the MGS mapped initialized
n_mapped_samples_vect <- c()
for(MGS in MGSs){
  if(paste(MGS,'_init',sep='')%in%names_avail){
    initial_genes <- trine_MGS_object[["i"]][[paste(MGS,'_init',sep='')]][1:100]
    
    if(paste(MGS,'_mean',sep='')%in%names_avail){
      mean_genes <- trine_MGS_object[["i"]][[paste(MGS,'_mean',sep='')]][1:100]
      
      #save the overlap
      overlap_table <- rbind(overlap_table,c(MGS,'init','mean',mean(initial_genes%in%mean_genes)))
      
      #if both of these exist and the best performance out exists, then save those
      if(paste(MGS,'_best',sep='')%in%names_avail){
        best_genes <- trine_MGS_object[["i"]][[paste(MGS,'_best',sep='')]][1:100]
        
        overlap_table <- rbind(overlap_table,c(MGS,'mean','best',mean(best_genes%in%mean_genes)))
        
        overlap_table <- rbind(overlap_table,c(MGS,'init','best',mean(best_genes%in%initial_genes)))
        
      }
    }
  }
  colsum <- colSums(utuseq_MGSgeneCountsML[[MGS]][1:n.genes, ])
  mapped <- colsum[colsum >= n_reads_minimum]
  
  n_mapped_samples_vect <- c(n_mapped_samples_vect,length(mapped))
}
#names(n_mapped_samples_vect) <- MGSs


#utuseq_MGSgeneCountsML = segataML
#load("/Users/trizac/Dropbox/CM/Artikel/R_Anders/coupon_simulations.Rdata")
load("/Users/trinezachariasen/Library/CloudStorage/OneDrive-DanmarksTekniskeUniversitet/CM/Artikel/R_Anders/coupon_simulations.Rdata")

if(file.exists(paste0(data_dir,'temp_MSE_table_filtered.Rdata')) & file.exists(paste0(data_dir,paste0('counts_list_table_filtered_',n_signature_genes_expected,'.Rdata')))){
  load(paste0(data_dir,'temp_MSE_table_filtered.Rdata'))
  load(paste0(data_dir,'counts_list_table_filtered_',n_signature_genes_expected,'.Rdata'))
}else{
  message('Recreating MSE table...')
  t0 <- Sys.time()
  #MSE table genertaionMGS
  MSE_table <- matrix(NA,ncol=3,nrow=0)
  colnames(MSE_table) <- c('MGS','stage','MSE')
  counts_list <- list()
  for(entry in names(trine_MGS_object[["i"]])){
    counts_list[[entry]] <- list()
    
    MGS <- paste(strsplit(entry,'_')[[1]][1],strsplit(entry,'_')[[1]][2],sep="_")
    
    #genes_of_interest <- names(present_genes[as.integer(trine_MGS_object[["i"]][[entry]][1:n.genes])])
    
    phase <- strsplit(entry,'_')[[1]][3]
    genes_of_interest <- outputs_parallel[[MGS]]$genes[phase]
    genes_of_interest <- genes_of_interest[[phase]]
    
    colsum <- colSums(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest,])
    mapped <- colsum[colsum >= 3]
    
    genes_r <- rbind(utuseq_MGSgeneCountsML[[MGS]][genes_of_interest, ])
    
    #then adjusts for gene length
    final.reads <- round(genes_r / ((present_genes[rownames(genes_r)]) * 10^-3))
    #Identifies number of genes unique signature genes per sample and the amount of counts assigned to those genes per sample
    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CCP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'fewer')) 
    }
    
    counts_list[[entry]][['pvals_analytical_fewer']] <- pvals
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CPP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'fewer')) 
    }
    counts_list[[entry]][['pvals_bootstrap_fewer']] <- pvals
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CPP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'fewer')) 
    }
    counts_list[[entry]][['pvals_fewer']] <- pvals
    
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CPP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'exact')) 
    }
    
    counts_list[[entry]][['pvals_analytical_exact']] <- pvals
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CPP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'exact')) 
    }
    counts_list[[entry]][['pvals_bootstrap_exact']] <- pvals
    
    pvals <- c()
    for(sample_name in names(countGenes)){
      pvals <-c(pvals,CPP_pval(k = countReads[sample_name],n = n_signature_genes_expected,d = countGenes[sample_name],mode = 'exact')) 
    }
    counts_list[[entry]][['pvals_exact']] <- pvals
    
    
    
    
    counts_list[[entry]][['countReads']] <- countReads
    #getting the number of identified signature genes for each sample
    counts_list[[entry]][['countGenes']] <- countGenes
    
    countPvals <- c()
    for(sample in names(countReads)){
      if(countReads[sample]>0){
        countPvals <- c(countPvals,CCP_pval(k = countReads[sample],n = as.numeric(n_signature_genes_expected),d = countGenes[sample],mode = 'fewer'))
      }else{
        countPvals <- c(countPvals,NA)
      }
    }
    counts_list[[entry]][['pval']] <- countPvals
    
    # calculating the predicted genecounts and then MSE
    pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
    
    # calculating the MSE
    mse <- mean((countGenes - pred)^2)
    
    MSE_table <- rbind(MSE_table,c(MGS,strsplit(entry,'_')[[1]][3],mse))
    message('Approximate time left:')
    print(((Sys.time()-t0)/
             match(entry,names(trine_MGS_object[["i"]])))*
            (length(trine_MGS_object[["i"]])-match(entry,names(trine_MGS_object[["i"]]))))
    
  }
  save(MSE_table,file = paste0(data_dir,'temp_MSE_table_filtered.Rdata'))
  save(counts_list,file = paste0(data_dir,paste0('counts_list_table_filtered_',n_signature_genes_expected,'.Rdata')))
  
}

#reformats
MSE_matrix <- matrix(NA,ncol=4,nrow=length(unique(MSE_table[,"MGS"])))
rownames(MSE_matrix) <- unique(MSE_table[,"MGS"])
colnames(MSE_matrix) <- c('init','mean','best','random')

for(row in seq(1,nrow(MSE_table))){
  MSE_matrix[MSE_table[row,"MGS"],MSE_table[row,"stage"]] <- as.numeric(MSE_table[row,"MSE"])
}

# Identifying the number of improved MSPs
improved <- length(MSE_matrix[,"best"][MSE_matrix[,"init"]-MSE_matrix[,"best"] > 0])

MGS_good_example <- "msp_54"
MGS_bad_example <- "msp_30"

#I find a lot of extremes when there are like 4 samples or whatever, which I don't think is enough to do jack, so I'm setting a limit. If an entry has fewer than minimum_sampels then fuck that
entries_with_n_or_more_samples <- c()
for(entry in names(counts_list)){
  if(sum(counts_list[[entry]][["countReads"]]>3)>minimum_sampels){
    entries_with_n_or_more_samples <- c(entries_with_n_or_more_samples,entry)
  }
}
#which MGSs have all 4 entries? init best and random
approved_MGSs <- table(unlist(strsplit(entries_with_n_or_more_samples, "([^_]+$)", perl=TRUE)))
names_of_MGS <- gsub('.{0,1}$', '', names(approved_MGSs)) #removing the last _
approved_MGSs <- sort(c(names_of_MGS[approved_MGSs==3],names_of_MGS[approved_MGSs==3],names_of_MGS[approved_MGSs==3]))
approved_entries <- paste0(approved_MGSs,c('_init','_mean','_best'))

MSE_matrix <- MSE_matrix[approved_MGSs,]

#trims overlap table
overlap_table <- overlap_table[overlap_table[,"MGS"]%in%approved_MGSs,]

#and the MSE table
MSE_table <- MSE_table[MSE_table[,"MGS"]%in%approved_MGSs,]

#table of PCCs between set1 and set2 signature gene set counts
PCC_table <- data.frame(MGS=character(),PCC=numeric(),stringsAsFactors = F)

for(MGS in names_of_MGS){
  PCC_table[nrow(PCC_table)+1,]=c(MGS,cor(counts_list[[paste0(MGS,'_',set1)]][['countReads']],
                                          counts_list[[paste0(MGS,'_',set2)]][['countReads']]))      
}
rownames(PCC_table) <- PCC_table[,"MGS"]


dims <- list(c('init','mean','best','random'),names_of_MGS,c('percent_failed','n_samples_failed','n_total_samples'))
failed_samples_array <- array(NA,dim = lengths(dims),dimnames = dims)
for(entry in approved_entries){
  MGS <- paste(strsplit(entry,'_')[[1]][1],strsplit(entry,'_')[[1]][2],sep="_")
  type <- strsplit(entry,'_')[[1]][3]
  failed_samples_array[type,MGS,"n_samples_failed"]=sum(counts_list[[entry]][["pval"]]<0.05,na.rm = T)
  failed_samples_array[type,MGS,"n_total_samples"]=sum(!is.na(counts_list[[entry]][["pval"]]),na.rm = T)
  failed_samples_array[type,MGS,"percent_failed"]=failed_samples_array[type,MGS,"n_samples_failed"]/failed_samples_array[type,MGS,"n_total_samples"]
}

message('Data generated')

################################################################################################
########################DATA GENERATED--------PLOTTING COMMENCES################################
################################################################################################

#plots correlation between number of mapped samples and overlap between initial and mean
plot_dat <- as.data.frame(
  cbind(n_mapped_samples_vect,
        as.numeric(overlap_table[overlap_table[,"set1"]=='init' & overlap_table[,"set2"]=='mean' ,"overlap"])))
colnames(plot_dat) <- c('n_sampled', 'overlap')
rownames(plot_dat) <- names_of_MGS

#only plots approved stuff

results <- cor.test(plot_dat[names_of_MGS,"n_sampled"],plot_dat[names_of_MGS,"overlap"],method='spearman',exact = NULL)
#plots correlation between number of mapped samples and overlap between initial and best
plot_dat <- as.data.frame(
  cbind(n_mapped_samples_vect,
        as.numeric(overlap_table[overlap_table[,"set1"]=='init' & overlap_table[,"set2"]=='best' ,"overlap"])))
colnames(plot_dat) <- c('n_sampled','overlap')
plot_dat <- plot_dat[approved_MGSs,]
results <- cor.test(plot_dat[,"n_sampled"],plot_dat[,"overlap"],method = 'spearman')


p <- ggplot(plot_dat,aes(n_sampled,overlap))+geom_point(aes(alpha=0.5))+
  xlab('Number of samples with MSP present')+
  ylab('Overlap between initial and final signature genes')+
  labs(title = 'Correlation between overlap and number of samples',subtitle = bquote("Spearman's rank correlation test p-value:" ~ .(signif(results[["p.value"]],3))~ rho~':'~.(signif(results[["estimate"]],2))))+
  theme(legend.position = 'none',plot.title = element_text(hjust=0.5,size=20),plot.subtitle = element_text(hjust=0.5,size=12))

#ggsave(p,filename = paste0(fig_dir,'initial_best_mapped_sample_dependnence.pdf'))


#plot that shows how correlated the overlap between two sets within an MGS and vs percentage-wise improvement in MSE
plotdat <- matrix(NA,ncol=5,nrow=0)
colnames(plotdat) <- c('MGS','overlap','MSE_initial','MSE_after','relative_MSE_improvement')
for(MGS in unique(approved_MGSs)){
  overlap <- overlap_table[overlap_table[,"MGS"]==MGS & overlap_table[,"set1"]==set1 & overlap_table[,"set2"]==set2,"overlap"]
  MSE_initial <- as.numeric(MSE_table[MSE_table[,"MGS"]==MGS & MSE_table[,"stage"]==set1,"MSE"])
  MSE_after <- as.numeric(MSE_table[MSE_table[,"MGS"]==MGS & MSE_table[,"stage"]==set2,"MSE"])
  relative_MSE_improvement <- MSE_after/MSE_initial
  
  plotdat <- rbind(plotdat,c(MGS,overlap,MSE_initial,MSE_after,relative_MSE_improvement))
}

plotdat <- as.data.frame(plotdat)
plotdat[,"overlap"] <- as.numeric(as.character(plotdat[,"overlap"]))
plotdat[,"MSE_initial"] <- as.numeric(as.character(plotdat[,"MSE_initial"]))
plotdat[,"MSE_after"] <- as.numeric(as.character(plotdat[,"MSE_after"]))
plotdat[,"relative_MSE_improvement"] <- as.numeric(as.character(plotdat[,"relative_MSE_improvement"]))

#saves some MGSs of special interest for later pointing out and sutff
MGSs_of_interest[["little_improvement_big_change"]] <- c(MGSs_of_interest[["little_improvement_big_change"]],as.character(plotdat[plotdat[,"relative_MSE_improvement"]<0.3 & plotdat[,"overlap"]<0.3,"MGS"]))
MGSs_of_interest[["big_improvement_little_change"]] <- c(MGSs_of_interest[["big_improvement_little_change"]],as.character(plotdat[plotdat[,"relative_MSE_improvement"]>0.75 & plotdat[,"overlap"]>0.75,"MGS"]))
MGSs_of_interest[["big_improvement_big_change"]] <- c(MGSs_of_interest[["big_improvement_little_change"]],as.character(plotdat[plotdat[,"relative_MSE_improvement"]>0.75 & plotdat[,"overlap"]<0.3,"MGS"]))

results <- cor.test(plotdat[, ][plotdat[, "overlap"] < 1, ]$overlap, plotdat[, ][plotdat[, "overlap"] < 1, ]$relative_MSE_improvement,method='spearman')

p <- ggplot(plotdat[,][plotdat[,'overlap']<1,],aes(overlap,relative_MSE_improvement))+geom_point(aes(alpha=0.5), show.legend = FALSE)+
  labs(x='Fraction of original signature genes retained',
       y=expression(paste("Relative MSE   ", frac(MSE [final], MSE [initial]))),
       title='Relative improvement in MSE as a function of initial signature genes retained',
       subtitle=bquote("Spearman's rank correlation test p-value:" ~ .(signif(results[["p.value"]],3))~ rho~':'~.(signif(results[["estimate"]],2))))+ theme_minimal() +
  theme(legend.position = 'none',plot.title = element_text(hjust=0.5, size=12),plot.subtitle = element_text(hjust=0.5)) + geom_density_2d() + geom_abline(intercept = 0, size=0.4)

p <- p + geom_label_repel(aes(label=MGS), show.legend = FALSE, box.padding = 0.35, point.padding=0.5) 

ggsave(p,filename = paste0(fig_dir,'overlap_vs_MSE_improvement.png'),height=7,width=7)

results <- cor.test(plotdat[,'overlap'],plotdat[,'MSE_initial'],method = 'spearman')
p <- ggplot(plotdat,aes(overlap,MSE_initial))+geom_point(aes(alpha=0.5))+
  labs(x='Fraction of original signature genes retained',y='Initial MSE',title='Initial MSE as a function of initial signature genes retained',subtitle = bquote("Spearman's rank correlation test p-value:" ~ .(signif(results[["p.value"]],3))~ rho~':'~.(signif(results[["estimate"]],2))))+
  theme(legend.position = 'none',plot.title = element_text(hjust=0.5,size=14),plot.subtitle = element_text(hjust=0.5)) + geom_density_2d()
ggsave(p,filename = paste0(fig_dir,'overlap_vs_init_MSE.png'), height=7,width=7)


#Plots how well the counts of signature genes before and after correlate with the number of exchanged sig. genes
temp <- overlap_table[overlap_table[,"set1"]==set1 & overlap_table[,"set2"]==set2,]
plotdat <- as.data.frame(cbind(PCC_table[temp[,1],"PCC"],temp[,4]))
colnames(plotdat) <- c('PCC','overlap')
plotdat[,"PCC"] <- as.numeric(as.character(plotdat[,"PCC"]))
plotdat[,"overlap"] <- as.numeric(as.character(plotdat[,"overlap"]))

cor_results <- cor.test(plotdat[,"PCC"],plotdat[,"overlap"],method = 'spearman')

#seems like some samples might largely change ID
p <- ggplot(plotdat,aes(overlap,PCC))+geom_point(aes(alpha=0.5))+scale_y_continuous(trans='log2')+
  labs(x='Fraction of original signature genes retained',
       y='PCC between original signature genes and refined signature genes',
       title='Correlation fraction of original signature genes retained and\n PCC between PCC between original and final signature gene set',
       subtitle=bquote("Spearman's rank correlation test p-value:" ~ .(signif(cor_results[["p.value"]],3))~ rho~':'~.(signif(cor_results[["estimate"]],3))))+
  theme(legend.position = 'none',plot.title = element_text(hjust=0.5,size=20),plot.subtitle = element_text(hjust=0.5,size=16)) 
ggsave(p,filename = paste0(fig_dir,'PCC_and_overlap.pdf'))



#saves interesting MGSs for later zooming in and stuff
MGSs_of_interest[["just_ignore_the_bad_samples"]] <- c(MGSs_of_interest[["just_ignore_the_bad_samples"]],MGSs[plotdat[,"PCC"]<0.8 & plotdat[,"overlap"]<0.3])


#Correlation between number of samples mapped with initial data set vs the change in sig. genes
#also gets data to see if a high degree of PCC also equates a high degree of rel MSE improvement
#also does a PCC correlation since I need all the same numbers anyway. Might as well add the PCC values and get two birds with one stone or whatever
samples_identified_vs_MSE <- data.frame(MGS=character(),relative_MSE_improvement=numeric(),relative_sample_change=numeric(),abs_sample_change=numeric(),PCC=numeric(),stringsAsFactors = F)
for(MGS in approved_MGSs){
  if(sum(MSE_table[,"MGS"]==MGS)>1){
    
    MSE_set1 <- as.numeric(MSE_table[MSE_table[,"MGS"]==MGS & MSE_table[,"stage"]==set1,"MSE"])
    MSE_set2 <- as.numeric(MSE_table[MSE_table[,"MGS"]==MGS & MSE_table[,"stage"]==set2,"MSE"])
    rel_MSE_improvement <- MSE_set2/MSE_set1
    
    n_samples_set1 <- sum(counts_list[[paste0(MGS,'_',set1)]][["countReads"]]>n_reads_minimum)
    n_samples_set2 <- sum(counts_list[[paste0(MGS,'_',set2)]][["countReads"]]>n_reads_minimum)
    
    abs_sample_change <- n_samples_set2-n_samples_set1
    rel_sample_change <- n_samples_set2/n_samples_set1
    
    samples_identified_vs_MSE[nrow(samples_identified_vs_MSE)+1,] <- c(MGS,
                                                                       rel_MSE_improvement,
                                                                       rel_sample_change,
                                                                       abs_sample_change,
                                                                       as.numeric(PCC_table[MGS,"PCC"]))
  }
  
}
old_col_names <- colnames(samples_identified_vs_MSE)

temp <- overlap_table[overlap_table[,"set2"]=='best' & overlap_table[,"set1"]=='init',]


#samples_identified_vs_MSE <- cbind(samples_identified_vs_MSE,as.numeric(temp[samples_identified_vs_MSE[,"MGS"],"overlap"]))
samples_identified_vs_MSE <- cbind(samples_identified_vs_MSE, as.numeric(temp[temp[,"MGS"] %in% samples_identified_vs_MSE[,"MGS"], 4]))
colnames(samples_identified_vs_MSE) <- c(old_col_names,'overlap')

samples_identified_vs_MSE[,"relative_MSE_improvement"] <- as.numeric(samples_identified_vs_MSE[,"relative_MSE_improvement"])
samples_identified_vs_MSE[,"relative_sample_change"] <- as.numeric(samples_identified_vs_MSE[,"relative_sample_change"])
samples_identified_vs_MSE[,"abs_sample_change"] <- as.numeric(samples_identified_vs_MSE[,"abs_sample_change"])
samples_identified_vs_MSE[,"PCC"] <- as.numeric(samples_identified_vs_MSE[,"PCC"])
#there appears to be a correlation between how many genes are excluded and how much better the fit becomes
results <- cor.test(samples_identified_vs_MSE[,"relative_sample_change"],samples_identified_vs_MSE[,"overlap"],method = 'spearman')

p <- ggplot(samples_identified_vs_MSE,aes(relative_sample_change,relative_MSE_improvement))+geom_point(aes(color=abs_sample_change),alpha=0.5)+
  scale_colour_gradientn("Absolute change in\nnumber of mapped\nsamples",colours = wesanderson::wes_palette("Zissou1", 1000, type = "continuous"))+
  theme(plot.title = element_text(hjust = 0.5,size=20),plot.subtitle = element_text(hjust=0.5))+
  labs(x= expression(paste("Relative number of mapped samples   ", frac(n ["samples final"], n ["samples initial"]))),
       y=expression(paste("Relative MSE   ", frac(MSE ["Final"], MSE ["Initial"]))),
       title = 'Change in MSE and change in mapped samples\nbefore and after signature gene refinement',
       subtitle=bquote("Spearman's rank correlation test p-value:" ~ .(signif(results[["p.value"]],3))~ rho~':'~.(signif(results[["estimate"]],2))),
       fill='Absolute change in mapped samples')+coord_flip()

ggsave(p,filename = paste0(fig_dir,'change_mapped_samples_and_MSE_change.pdf'))


p_box <- qplot(x = 1,y=samples_identified_vs_MSE[,"relative_sample_change"],geom='violin')+geom_boxplot(width=0.05,outlier.size = 0.7,outlier.alpha = 0.4)+theme(axis.text.x = element_blank())+
  labs(y=expression(paste("Relative number of mapped samples   ", frac(n ["samples final"], n ["samples initial"]))),title='Distribution of relative number of mapped samples')



title <- ggdraw() + draw_label("Change in number of mapped samples", fontface='bold')

saveme <-cowplot::plot_grid(p_box+theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.title=element_blank(),plot.margin = unit(c(0,0.2,0.4,0),'inches') ),
                            p+theme(axis.title.y = element_blank(),axis.text.y = element_text(),axis.text.y.left = element_blank(),axis.ticks.y.left=element_blank(),plot.margin=unit(c(0,0,0,0.2),'inches'),
                                    plot.title = element_blank(),plot.subtitle = element_blank()),
                            nrow=1,
                            labels=c('A','B'),
                            rel_widths = c(1,3))

saveme <- cowplot::plot_grid(title,saveme,ncol=1,rel_heights=c(0.1,1))

ggsave(saveme,filename = paste0(fig_dir,'test.pdf'),height=7,width=7)

#saves MGSs of interest for later focusing
MGSs_of_interest[["just_ignore_the_bad_samples"]] <- c(MGSs_of_interest[["just_ignore_the_bad_samples"]],samples_identified_vs_MSE[samples_identified_vs_MSE[,"relative_MSE_improvement"]<(-0.75) & samples_identified_vs_MSE[,"abs_sample_change"]<(-200),"MGS"])

#to which degree is the improvement in MSE just a matter of changing the identity of the thing measured? Not so much I think
#to some extend sure and I think we need to consider imposing a cut-off value for the degree of persistence of identity.
#not a problem for most part though

results <- cor.test(samples_identified_vs_MSE[,"PCC"],samples_identified_vs_MSE[,"relative_MSE_improvement"],method = 'spearman')
p <- ggplot(samples_identified_vs_MSE,aes(PCC,relative_MSE_improvement))+geom_point(alpha=0.5)+
  labs(x='PCC between old and new set of genes',y=expression(paste("Relative MSE   ", frac(MSE [final], MSE [initial]))),
       title='Correlation of initial and final signature genes',
       subtitle=bquote("Spearman's rank correlation test p-value:" ~ .(signif(results[["p.value"]],3))~ rho~':'~.(signif(results[["estimate"]],2))))+theme(plot.title = element_text(size = 20,hjust=0.5),plot.subtitle=element_text(hjust=0.5))

ggsave(p,filename = paste0(fig_dir,'PCC_and_rel_MSE_improvement.pdf'))

MGSs_of_interest[["just_ignore_the_bad_samples"]] <- c(MGSs_of_interest[["just_ignore_the_bad_samples"]],samples_identified_vs_MSE[samples_identified_vs_MSE[,"relative_MSE_improvement"]<(-0.75) & samples_identified_vs_MSE[,"PCC"]<0.8,"MGS"])

#MSEs before and after plot
init_values <- as.numeric(MSE_table[MSE_table[,"stage"]=='init',"MSE"])
#step1_values <- as.numeric(MSE_table[MSE_table[,"stage"]=='mean',"MSE"])
step2_values <- as.numeric(MSE_table[MSE_table[,"stage"]=='best',"MSE"])

plotdat <- as.data.frame(
  rbind(cbind(init_values,rep('Initial',length(init_values))),
#        cbind(step1_values,rep('Step 1',length(step1_values))),
        cbind(step2_values,rep('Final',length(step2_values)))
  )
)
colnames(plotdat) <- c('value','type')
plotdat[,"value"] <- as.numeric(as.character(plotdat[,"value"]))


#library(ggpubr) 
p <- ggplot(plotdat,aes(type,value))+geom_boxplot()+scale_y_continuous(trans = 'log',breaks = scales::log_breaks())+stat_compare_means(paired = TRUE, position=position_jitter(width=0, height=0.1))+
  labs(x='Refinement stage',y='MSE',title='MSEs across steps')+theme(plot.title = element_text(hjust=0.5,size=20))
p
ggsave(p,file=paste0(fig_dir,'MSE_before_after_boxplot.pdf'))


#alternate dot-plot
plotdat <- as.data.frame(cbind(MSE_table[MSE_table[,"stage"]=='init',c("MGS","MSE")],MSE_table[MSE_table[,"stage"]=='best',"MSE"]))

colnames(plotdat) <- c('MGS','MSE_init','MSE_best')
plotdat[,"MSE_init"] <- as.numeric(as.character(plotdat[,"MSE_init"]))
plotdat[,"MSE_best"] <- as.numeric(as.character(plotdat[,"MSE_best"]))

range <- c(min(plotdat[,c("MSE_init","MSE_best")]),max(plotdat[,c("MSE_init","MSE_best")]))
results <- wilcox.test(plotdat[,"MSE_init"],plotdat[,"MSE_best"])
p <- ggplot(plotdat,aes(MSE_init,MSE_best))+geom_point(alpha=0.5)+
  #log scale
  scale_y_continuous(trans = 'log',breaks = scales::log_breaks(),limits = range)+
  scale_x_continuous(trans = 'log',breaks = scales::log_breaks(),limits = range)+
  geom_abline(intercept=0,slope = 1,color=cmcol1[1])+
  labs(x='Initial MSE',y='Final MSE',title='MSEs before and after two-step refinement',
       subtitle=paste0('Wilcoxon rank test initial vs final p-value: ',format(signif(results[["p.value"]],3), scientific=TRUE)))+theme(plot.title = element_text(hjust=0.5,size=20),plot.subtitle = element_text(hjust=0.5)) + theme_minimal()
ggsave(p,filename = paste0(fig_dir,'MSE_before_after.png'))
MGSs_of_interest[["big_change_MSE"]] <- c(MGSs_of_interest[["big_change_MSE"]],as.character(plotdat[plotdat[,"MSE_best"]-plotdat[,"MSE_init"]==min(plotdat[,"MSE_best"]-plotdat[,"MSE_init"]),"MGS"]))

# The number of MSPs which has a lover MSE after refinement
plotdat[,"MSE_best"][(plotdat[,"MSE_init"]-plotdat[,"MSE_best"])!=0]
# Total number of MSPs
length(plotdat[,"MSE_best"])





#a plot that illustrates the pdf of interest
#plots for all number of samples
#and all d's
if(plot_general_illustration==T){
  
  
  pvals_boot <- c()
  pvals_analytical <- c()
  
  plotdat <- matrix(NA,ncol=3,nrow=0)
  for(k in seq(1,1900)){
    print(k)
    for(d in seq(1,n.genes)){
      plotdat <- rbind(plotdat,c(k,d,coupon_sim_pval_analytical_if_possible(n_reads = k,n_signature_genes_expected = '100',observed_signature_genes = d,mode = 'exact')))
      pvals_boot <- c(pvals_boot,coupon_sim_pval_just_bootstrap(n_reads = k,n_signature_genes_expected = '100',observed_signature_genes = d,mode = 'exact'))
      pvals_analytical <- c(pvals_analytical,coupon_sim_pval_just_analytical(n_reads = k,n_signature_genes_expected = '100',observed_signature_genes = d,mode = 'exact'))
    }
  }
  
  
  #plots bootstrapped p-vals vs analytical p-vals
  names <- c()
  for(entry in names(counts_list)){
    names <- c(names,paste0(entry,'-',names(counts_list[[entry]][["countReads"]])))
    pvals_boot <- c(pvals_boot,counts_list[[entry]][["pvals_bootstrap_exact"]])
    pvals_analytical <- c(pvals_analytical,counts_list[[entry]][["pvals_analytical_exact"]])
  }
  plotdat <- data.frame(cbind(names[!is.na(pvals_boot) & !is.na(pvals_analytical)],
                              pvals_boot[!is.na(pvals_boot) & !is.na(pvals_analytical)],
                              pvals_analytical[!is.na(pvals_boot) & !is.na(pvals_analytical)]))
  
  colnames(plotdat) <- c('sample','boot','analytical')
  plotdat[,"boot"] <- as.numeric(as.character(plotdat[,"boot"]))
  plotdat[,"analytical"] <- as.numeric(as.character(plotdat[,"analytical"]))
  PCC <- cor(plotdat[,"boot"],plotdat[,"analytical"],method='pearson')
  
  
  p <- ggplot(plotdat,aes(analytical,boot))+geom_point(alpha=0.01)+xlab('True value')+ylab('Bootstrapped value')+labs(title='Comparison of bootstrapped value and true values for all tests',subtitle=paste0('PCC: ',as.character(round(PCC,3))))+
    geom_abline(intercept = 0,slope = 1)+theme(plot.title=element_text(hjust=0.5),
                                               plot.subtitle = element_text(hjust=0.5))
  ggsave(p,filename = paste0(fig_dir,'pvalues_boot_vs_true.pdf'))
  
  #plots an illustrative eplot of the PDFs
  colnames(plotdat) <- c('Reads','Genes','p')
  plotdat <- as.data.frame(plotdat)
  p <- ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
              data = plotdat, aes(y = Genes, x = Reads)) + 
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes), size=1, alpha=0.2,colour="#b8dde3") + #the expected distribution
    geom_point(size = 0.7, aes(color = p)) +
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected')

  p <- (p + scale_x_log10())
  
  p <- p  + theme_minimal() + theme(legend.position="bottom")+scale_colour_gradientn("Chance of complete\n signature gene set",colours = wesanderson::wes_palette("Darjeeling1", 1000, type = "continuous"))
  
 # ggsave(p,filename = paste0(fig_dir,'exact_chance_100.pdf'))
  
  #aaand for fewer. Because I think people will like that more
  plotdat <- matrix(NA,ncol=3,nrow=0)
  for(k in seq(1,1900)){
    #print(k)
    for(d in seq(1,n.genes)){
      plotdat <- rbind(plotdat,c(k,d,coupon_sim_pval_analytical_if_possible(n_reads = k,n_signature_genes_expected = '100',observed_signature_genes = d,mode = 'fewer')))
    }
  }
  
  
  colnames(plotdat) <- c('Reads','Genes','p')
  
  plotdat[plotdat[,"p"]<0.05,"p"] <- NA
  plotdat <- as.data.frame(plotdat)
  p <- ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
              data = plotdat, aes(y = Genes, x = Reads)) + 
    geom_point(size = 0.7, aes(color = p)) +
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes), size=1, alpha=.6,colour="#000000") + #t he expected distribution
    
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected')
  
  p <- (p + scale_x_log10())
  p <- p  + theme_minimal() + theme(legend.position="bottom")+scale_colour_gradientn("Chance of complete\n signature gene set",colours = wesanderson::wes_palette("Darjeeling1", 1000, type = "continuous"),na.value = 'black')
  
  
  p
  #ggsave(p,filename = paste0(fig_dir,'fewer_chance_100.pdf'))
  
  
  
  
  weird_sampels_names <- as.character(plotdat[plotdat[,"boot"]>0.65 & plotdat[,"boot"]<0.7,"sample"])
  
  weird_samples_n <-c()
  weird_samples <- list()
  for(i in weird_sampels_names){
    entry <- strsplit(i,'[-]')[[1]][1]
    sample <- strsplit(i,'[-]')[[1]][2]
    weird_samples[[i]] <- c(counts_list[[entry]][['countGenes']][sample],counts_list[[entry]][['countReads']][sample],counts_list[[entry]][["pvals_bootstrap"]][sample])
    weird_samples_n <- c(weird_samples_n,counts_list[[entry]][['countReads']][sample])
  }
}

#plots that plot that lets us adjust the n-level
#gathers all the plotdat first
plotdat <- data.frame()
backup_sig_genes_expected <- n_signature_genes_expected

MGS_good_example <- "msp_54"
MGS_bad_example <- "msp_11"

for(n_signature_genes_expected in c(85,90,95,100)){
  #gathers plotdat for an MGS with a good distribution at a certain threshold
  plotdat_temp <- pois.plotdata.with.probabilities(Genes = counts_list[[paste0(MGS_good_example,'_init')]][['countGenes']],
                                                   Reads = counts_list[[paste0(MGS_good_example,'_init')]][['countReads']],
                                                   mse =  as.numeric(MSE_table[MSE_table[,"MGS"]==MGS_good_example & MSE_table[,"stage"]=="init","MSE"]),
                                                   id = MGS_good_example,HGMGS = HGMGS,
                                                   mode = 'exact',n_signature_genes_expected = n_signature_genes_expected) #fewer
  plotdat_temp <- cbind(plotdat_temp,rep(n_signature_genes_expected,nrow(plotdat_temp)),rep(MGS_good_example,nrow(plotdat_temp)))
  colnames(plotdat_temp) <- c('Reads','Genes','pval','name','n_lvl','MGS')
  plotdat <- rbind(plotdat,plotdat_temp)
  
  #gathers plotdat for an MGS with a bad distribution at a certain threshold
  plotdat_temp <- pois.plotdata.with.probabilities(Genes = counts_list[[paste0(MGS_bad_example,'_init')]][['countGenes']],
                                                   Reads = counts_list[[paste0(MGS_bad_example,'_init')]][['countReads']],
                                                   mse = as.numeric(MSE_table[MSE_table[,"MGS"]==MGS_bad_example & MSE_table[,"stage"]=="init","MSE"]),id = MGS_bad_example,HGMGS = HGMGS,
                                                   mode = 'exact',n_signature_genes_expected = n_signature_genes_expected)
  plotdat_temp <- cbind(plotdat_temp,rep(n_signature_genes_expected,nrow(plotdat_temp)),rep(MGS_bad_example,nrow(plotdat_temp)))
  colnames(plotdat_temp) <- c('Reads','Genes','pval','name','n_lvl','MGS')
  plotdat <- rbind(plotdat,plotdat_temp)
}
plotdat[,"n_lvl"] <- factor(plotdat[,"n_lvl"],rev(c(85,90,95,100)))

n_signature_genes_expected<- backup_sig_genes_expected


p = ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
           data = plotdat, aes(text = name, y = Genes, x = Reads)) + 
  stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes), size=1, alpha=0.2,colour="#b8dde3") + #t he expected distribution
  geom_point(size = 0.3, aes(color = pval)) +
  xlab('Reads mapped (log)') +
  ylab('Signature genes detected') +
  
  ggtitle(paste0('Performances of '))#,taxonomy$MainTax[taxonomy$Cluster==MGS_bad_example],' and\n',taxonomy$MainTax[taxonomy$Cluster==MGS_good_example],' across thresholds'))

  #ggtitle(paste0('Performances of ',HGMGS[["Tax"]][MGS_bad_example,"MainTax"],' and\n',HGMGS[["Tax"]][MGS_good_example,"MainTax"],' across thresholds'))

p <- (p + scale_x_log10())
p <- p  + theme(legend.position="bottom",plot.title=element_text(hjust=0.5,size=20))+scale_colour_gradientn("Chance of complete\n signature gene set",colours = wesanderson::wes_palette("Darjeeling1", 1000, type = "continuous"))
p <- p+facet_grid(n_lvl~MGS)
p
#ggsave(p,filename = paste0(fig_dir,'thresholds_test_exact.pdf'),width=10,height=12)



#plots every effing poisson distribution. Takes a while!
if(plot_all_poissons ==T){
  print('plots every single possion distribution')
  i=0
  #plots all MGSs before/during/after refinement in a pdf in each three 
  pdf(paste0(fig_dir,'ALL_DISTRIBUTIONS_exact.pdf'))
  #pdf(paste0(fig_dir,'Clusters_no_improvement_DISTRIBUTIONS.pdf'))
  #pdf(paste0(fig_dir,'Clusters_large_decrease_mapped_samples_DISTRIBUTIONS.pdf'))
  
  #Clusters_no_improvement <- c("Cluster33","Cluster40","Cluster7124","Cluster343","Cluster433","Cluster724","Cluster179","Cluster635","Cluster222","Cluster403","Cluster408","Cluster347","Cluster895","Cluster314","Cluster10707","Cluster112","Cluster2143","Cluster83","Cluster1112","Cluster237","Cluster37","Cluster14","Cluster658","Cluster1203","Cluster244","Cluster404")
  #min_20_relative_mapped_samples <- samples_identified_vs_MSE[order(samples_identified_vs_MSE$relative_sample_change),"MGS"][1:20]
  
  for(MGS in approved_MGSs){
 # for(MGS in Clusters_no_improvement){
    i=i+1
    plotdat <- data.frame()
    for(set in c('init','best')){
      if(set=='init'){
        print_type <- 'Initial'
#      }else if(set=='mean'){
 #       print_type <- 'Step1'
      }else if(set=='best'){
        print_type <- 'Final'
      }
      
      
      
      plotdat_temp <- pois.plotdata.with.probabilities(Genes = counts_list[[paste0(MGS,'_',set)]][['countGenes']],
                                                       Reads = counts_list[[paste0(MGS,'_',set)]][['countReads']],
                                                       mse =  as.numeric(MSE_table[MSE_table[,"MGS"]==MGS & MSE_table[,"stage"]==set,"MSE"]),
                                                       id = MGS,
                                                       HGMGS = HGMGS,
                                                       mode = 'exact',
                                                       n_signature_genes_expected = n_signature_genes_expected)
      plotdat_temp <- cbind(plotdat_temp,rep(print_type,nrow(plotdat_temp)))
      colnames(plotdat_temp) <- c('Reads','Genes','pval','name','type')
      plotdat <- rbind(plotdat,plotdat_temp)
    }
    
    tax <- tax_df$Taxonomy[tax_df$MSP == MGS]
    
    p = ggplot(log = "x", xlim = c(1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
               data = plotdat, aes(text = name, y = Genes, x = Reads)) + 
      stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes), size=1, alpha=0.2,colour="#b8dde3") + #t he expected distribution
      geom_point(size = 0.3, aes(color = pval)) +
      xlab('Reads mapped (log)') +
      ylab('Signature genes detected') +
      ggtitle(paste0('Performances of ',MGS, ': ', tax, ''))

    p <- (p + scale_x_log10(limits=c(1,100000)))
    p <- p  + theme(legend.position="bottom")+scale_colour_gradientn("Chance of complete\n signature gene set",colours = wesanderson::wes_palette("Darjeeling1", 1000, type = "continuous"))
    p <- p+facet_grid(~factor(type, levels=c('Initial', "Step1", "Final")))
    p <- p+labs(subtitle=paste0('Initial MSE: ',round(MSE_matrix[MGS,'init'],2),' Refined MSE: ',round(MSE_matrix[MGS,'best'],2)))
    
    print(p)
    accepted_samples <- table(na.omit(plotdat)[na.omit(plotdat)$pval>0.05,"type"]) # the number of samples, that are accepted for the MGS for the initial and final SGs
    rejected_samples <- table(na.omit(plotdat)[na.omit(plotdat)$pval<0.05,"type"]) # the number of samples, that are rejected for the MGS
    print(c(accepted_samples, rejected_samples))
  }
  dev.off()
}

 #frequency distribution


plotdat <- matrix(NA,nrow=0,ncol=3)
colnames(plotdat) <- c('index','frequency','type')
#plots hetereogeneous frequencies of sets for MGS0010
MGS_of_interest <- 'Cluster10'
for(entry in c('init','mean','best')){
  signature_genes <- trine_MGS_object[["i"]][[paste0(MGS_of_interest,'_',entry)]][1:100]
  signature_genes <- names(present_genes[signature_genes])
  frequency <- rowSums(utuseq_MGSgeneCountsML[[MGS_of_interest]][signature_genes,]>0)/ncol(utuseq_MGSgeneCountsML[[MGS_of_interest]][signature_genes,])
  frequency <- sort(frequency,decreasing = T)
  names(frequency) <- seq(1,length(frequency))
  
  plotdat_temp <- cbind(names(frequency),frequency,rep(entry,length(frequency)))
  colnames(plotdat_temp) <- c('index','frequency','type')
  plotdat <- rbind(plotdat,plotdat_temp)
}

plotdat[,"type"] <- gsub('best','Final',gsub('mean','Step 1',gsub('init','Initial',plotdat[,"type"])))

plotdat <- as.data.frame(plotdat)

plotdat[,'type'] <-factor(as.character(plotdat[,"type"]),c('Initial','Step 1','Final')) 

plotdat[['index']] <- factor(as.numeric(as.character(plotdat[,"index"])),
                             levels = sort(unique(as.numeric(as.character(plotdat[,"index"]))),decreasing = F))

#plotdat[["index"]] <- as.numeric(as.character(plotdat[["index"]]))
plotdat[["frequency"]] <- as.numeric(as.character(plotdat[["frequency"]]))
title_str <- paste0('Frequencies of signature genes detection for ', MGS_of_interest)

p <- ggplot(data = plotdat,aes(x = index,y = frequency))+geom_bar(stat='identity')+facet_grid(~type)+theme_minimal()+
  labs(title = title_str ,x='Gene index',y='Share of samples detected in')+
  theme(axis.text.x=element_blank(),plot.title = element_text(hjust=0.5,size=14))

ggsave(p,filename = paste(fig_dir, MGS_of_interest, "_frequencies.pdf", sep=""))


#fraction of samples rejected before and after
rejected_percent <- MSE_matrix
rejected_percent[,] <- NA

for(MGS in approved_MGSs){
  for(set in c('init','mean','best','random')){
    entry <- paste0(MGS,'_',set)
    rejected_percent[MGS,set] <- mean(percent_rejected <- counts_list[[entry]][["pval"]]<0.05,na.rm = T)*100
  }
}



#alternate dot-plot
plotdat <- as.data.frame(cbind(rejected_percent[,"init"],rejected_percent[,"best"],rownames(rejected_percent)))

colnames(plotdat) <- c('rejected_init','rejected_final','MGS')
plotdat[,"rejected_init"] <- as.numeric(as.character(plotdat[,"rejected_init"]))
plotdat[,"rejected_final"] <- as.numeric(as.character(plotdat[,"rejected_final"]))


range <- c(1,100)
results <- wilcox.test(plotdat[,"rejected_init"],plotdat[,"rejected_final"])
p <- ggplot(plotdat,aes(rejected_init,rejected_final))+geom_point(alpha=0.5)+ geom_density_2d()+
  #log scale
  scale_y_continuous(trans = 'log',breaks = scales::log_breaks(),limits = range)+
  scale_x_continuous(trans = 'log',breaks = scales::log_breaks(),limits = range)+
  geom_abline(intercept=0,slope = 1,color=cmcol1[1])+
  labs(x='Percent rejected for initial signature genes',
       y='Percent rejected for refined signature genes',
       title='Percent of rejected samples per MGS',
       subtitle=bquote("Spearman's rank rank sum test p-value:" ~ .(signif(results[["p.value"]],3))))+theme(plot.title = element_text(hjust=0.5,size=20),plot.subtitle=element_text(hjust=0.5))

ggsave(p,filename = paste0(fig_dir,'Percent_rejected_before_after.pdf'))
#rejected samples as function of MSE


#alternate dot-plot
plotdat <- as.data.frame(rbind(cbind(rejected_percent[,"best"],rep('Final',nrow(rejected_percent))),
                               cbind(rejected_percent[,"init"],rep('Initial',nrow(rejected_percent)))))


colnames(plotdat) <- c('value','type')
plotdat[,"type"] <- factor(as.character(plotdat[,"type"]),c('Initial','Final'))
plotdat[,"value"] <- as.numeric(as.character(plotdat[,"value"]))

result <- wilcox.test(plotdat[plotdat[,"type"]=='Initial',"value"],plotdat[plotdat[,"type"]=='Final',"value"])
p <- ggplot(plotdat,aes(type,value))+geom_boxplot()+scale_y_continuous(trans = 'log',breaks = scales::log_breaks())+
  labs(x='Refinement stage',y='Percent rejected',title='Percent rejected before and after refinement',subtitle=paste0('Wilcoxon test p-value: ',as.character(signif(result[['p.value']],3))))+
  theme(plot.title = element_text(hjust=0.5,size=20),plot.subtitle = element_text(hjust=0.5,size=15))
ggsave(p,file=paste0(fig_dir,'rejected_percent_before_after_boxplot.pdf'))







############stats for the article############
temp <- unlist(counts_list)
print(paste0('the mean of samples that pass <0.05 at initial'))
pval_init <- temp[grepl('_init',names(temp)) & grepl('pval',names(temp))]
print(mean(pval_init<0.05,na.rm = T))

mean(temp[temp[grepl('_init',names(temp)) & grepl('pval',names(temp))]>0.05 & grepl('pval',names(temp))]>0.05,na.rm=T)
print(paste0('the mean of samples that pass <0.05 at mean'))
pval_mean <- temp[grepl('_mean',names(temp)) & grepl('pval',names(temp))]
print(mean(pval_mean<0.05,na.rm = T))


print(paste0('the best of samples that pass <0.05 at best'))
pval_best <- temp[grepl('_best',names(temp)) & grepl('pval',names(temp))]
print(mean(pval_best<0.05,na.rm = T))


print("number of total attempted MGSs")
print(length(tax_df$MSP))
print("number of approved MGSs")
print(length(entries_with_n_or_more_samples)/3)
print('discarded MGSs')
print(length(tax_df$MSP)-length(entries_with_n_or_more_samples)/3)


###################
# Getting the number of mapped samples and number of Genes within the Clusters, which do not improve during refinement

# 
# Clusters_no_improvement <- c("Cluster33","Cluster40","Cluster7124","Cluster343","Cluster433","Cluster724","Cluster179","Cluster635","Cluster222","Cluster403","Cluster408","Cluster347","Cluster895","Cluster314","Cluster10707","Cluster112","Cluster2143","Cluster83","Cluster1112","Cluster237","Cluster37","Cluster14","Cluster658","Cluster1203","Cluster244","Cluster404")
# 
# avg_mapped <- 0
# avg_genes <- 0
# counter <- 0
# 
# for (MGS in names(Clusterlist)){  #Clusters_no_improvement){
#   genes <- dim(Clusterlist[[MGS]])[1]
#   mapped_samples <- length(counts_list[[paste0(MGS,'_',"init")]][['countGenes']][counts_list[[paste0(MGS,'_',"init")]][['countGenes']]>0])
#  # print(c(MGS, genes, mapped_samples))
#   counter <- counter + 1
#   
#   avg_mapped <- avg_mapped + mapped_samples
#   avg_genes <- avg_genes + genes
#   
# }
# 
# avg_mapped/counter # The average number of mapped samples per cluster
# avg_genes/counter  # The average number of genes per cluster
# 
# #### Examining the genes, that we suspect are responsible for the crossmapping 
# #### / the pattern we see in the 20 samples, with the largest decrease in number of mapped samples
# MGS_object_dir <- paste(data_dir,'VAMB_refined.RDS',sep='')
# outputs_parallel <- readRDS(file=MGS_object_dir)
# 
# discarded_list <- c()
# for (cluster in min_20_relative_mapped_samples){
#   init_genes <- outputs_parallel[[cluster]]$genes$init
#   final_genes <- outputs_parallel[[cluster]]$genes$best
#   
#   discarded_init <- init_genes[!init_genes %in% final_genes]
#   
#   discarded_list <- c(discarded_list, discarded_init)
# }
# 
# # saving the discarded genes to a file
# write.table(discarded_list, paste(work_dir, "discarded_genes.txt", sep=""), row.names=F, col.names=FALSE,quote=FALSE)
# 
# 











