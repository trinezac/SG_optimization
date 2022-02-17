# loading functions

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
  
  replace.genes <- sort(abs(gene.prob) * 0.975 + gene.pos * 0.025, decreasing = FALSE)[0:n.replace]
  
  return("replace" = replace.genes)
}

rank.mat <- function(data, id, gene.names, step, phase, threshold, t.start, t.end, n.replace = FALSE, mapped) {
  #  Identify suitable signature genes by fitting NB model and ranking the genes
  #
  # Args:
  #   data: The object containing the biological entities and genes
  #   gene.names: A list of genes used to train the model
  #   step: Indicating whether the suitable genes are found based on "median" or "quantile"
  #   phase: The phase of the modelling. Currently there are two phases: 
  #     "initial": fit NB model to the original signature genes
  #     "rotation": fit NB model to refined signature genes 
  #   threshold: Mean value of rank for keeping genes
  #   t.start: The start of the genes in the read count matrix, default: 0 for initial, 100 for rotation
  #   t.end: The end of the genes in the read count matrix
  #   n.replace: The number of genes to replace, int
  #   mapped: samples we are confident the entity is present in (e.g. have more than n.mapped.minimum reads mapping)
  #
  # Returns:
  #   rank.mat: Matrix containing the ranks of the gene.names in the samples
  #   good.genes: List of the gene names of the genes that we keep in the signature gene set, for rotation it is the list of genes of length n.replace
  #   countReads: double with sample name and corresponding total read count
  #   countGenes: double with sample name and corresponding detected number of signature genes
  #   mse: The MSE of the fit to the expected distribution (s=(1-((G-1)/G)^N)*G)
  
  rank.matrix <- matrix(NA, nrow=(t.end-t.start), ncol=length(mapped))
  
  reads <- round(data[[id]][gene.names, names(mapped)] / 
                   (present_genes[rownames(data[[id]][gene.names, ])] * 10^-3))
  
  # Fitting the negative binomial distribution on all samples, saving the model in the gf matrix
  gf.nbin <- apply(matrix(1:length(mapped), nrow = length(mapped), ncol = 1), 1, function(x) (goodfit(reads[, x], type = "nbinom", method = "ML", par = list(size = mean(reads[,x])))))
  
  
  # filling in the rank.matrix with the rank of all genes within each sample
  if (phase == "initial") {
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[reads[,c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
    
  } else if (phase == "rotation") {
    new.reads <- matrix(NA, nrow = (t.end-t.start), ncol = length(mapped))
    new.reads <- round(data[[id]][(t.start + 1):t.end, 1:length(names(mapped))] / 
                         (present_genes[rownames(data[[id]][(t.start + 1):t.end, ])] * 10^-3))
    
    # if a new gene have more reads than found in the good genes, then set the readcount to the max observed in the good genes
    for (i in 1:length(mapped)) new.reads[, i][new.reads[, i] > max(gf.nbin[[i]]$count)] = max(gf.nbin[[i]]$count)
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[new.reads[, c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
  }
  
  rownames(rank.matrix) <- names(data[[id]][(t.start + 1):t.end, 1]) # genes as rownames
  colnames(rank.matrix) <- names(mapped) # samples as colnames
  
  #######################################################################
  # Adding step to weigh samples with more reads higher
  ######################################################################
  ori_mean = mean(rowMeans(rank.matrix)) #50.5 
  min = min(rowMeans(rank.matrix))- 5
  max = max(rowMeans(rank.matrix)) + 5
  rank.matrix <- t(t(rank.matrix) * log10(colSums(reads))) 
  rank.matrix <- replace(rank.matrix, rank.matrix == -Inf, 0)
  
  scale = ori_mean/mean(rowMeans(rank.matrix))
  rank.matrix <- replace(rank.matrix*scale, rank.matrix*scale == 'NaN', 0)

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
    
    genes_r <- data[[id]][c(names(good.genes),gene.names), ]
    
    
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))

    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
  } else if (phase == "rotation"){
    if (step == "mean"){
      gene.performance <- rowMeans(rank.matrix)
    } else if(step=="quantile"){
      gene.performance <- rowQuantiles(rank.matrix, probs = seq(from = 0, to = 1, by = 0.05))[,20] #95percentile
    }
    
    names(gene.performance) = names(data[[id]][(t.start + 1):t.end, 1]) 
    high.rank <- gene.performance 
    high.pos <- seq((t.start + 1), t.end)
    
    # ensuring that the new genes are not already part of the signature gene set
    high.rank <- high.rank[! names(gene.performance) %in% gene.names]
    high.pos <- high.pos[! names(gene.performance) %in% gene.names]
    
    # identifying the best replacement genes
    good.genes <- gene.selection(n.replace, high.rank, high.pos)
    
    genes_r <- data[[id]][c(names(good.genes),gene.names), ]
    
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))

    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
  }
  
  # calculating the predicted gene counts
  pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
  
  # calculating the MSE
  mse <- mean((countGenes - pred)^2)
  
  return(list("good.genes" = good.genes, "countReads" = countReads, "countGenes" = countGenes, "mse" = mse,'countReads'=countReads,'countGenes'=countGenes))
}

run_one_id <- function(id){
  # For each MSP go through the 100 inital SGs, identifying genes to replace
  # New genes are tested and inserted to the SG set if they entail an improved performance 
  # The final SG set is returned along with it's MSE compared to the expected distribution
  
  # Args:
  #   id: The ID of the biological entity
  # Returns:
  #   outputlist: a list containing the: 
  #        'id': id
  #        'MSE': a list of the MSE of the initial, mean and final SG sets
  #        'genes': a list of lists, containing the 100 genes from the initial, mean and final sets
  
  outputlist <- list()
  outputlist[['id']] <- id  # the id of the biological entity
  outputlist[['MSE']] <- list()
  outputlist[['genes']] <- list()
  message(id)
  
  #check that more than n genes represent in the id
  if (nrow(MSPlist[[id]])<=n.genes){
    print(c("Not enough genes in the MSP The number is", nrow(MSPlist[[id]])))
    return()
  }
  
  #colsum is the amount mapped genes in any given sample
  colsum <- colSums(MSPlist[[id]][1:n.genes, ])
  # extracting only the samples which contains above 3 mapped reads -- if there are 3 reads we believe this is a true detection
  mapped <- colsum[colsum >= n.mapped.minimum]
  
  n_samples <- setNames(c(n_samples, sum(colsum>=n.mapped.minimum)), c(names(n_samples), id))
  
  #if less than 3 samples are mapped to the MGS, the MGS should be skipped
  if (length(mapped) < 3){
    return()}
  
  genes <- names(MSPlist[[id]][1:n.genes, 1])
  best.genes.step.zero <- genes
  
  ####################################################### The mean gene refinement #########################################################
  # loop over different threshold values for mean rank
  best.genes <- c()
  mse <- c()
  
  #fit without cutting anything out at all
  init.fit <- rank.mat(MSPlist, id = id,
                       gene.names = genes,
                       step = "mean",
                       phase = "initial",
                       threshold = 999,
                       t.start = 0,
                       t.end = n.genes,
                       mapped = mapped)
  
  #special case:
  #cant fix perfect. just stop here if you come to that 
  if(init.fit[['mse']]==0){
    outputlist[['MSE']][['initial']] <- 0
    outputlist[['MSE']][['mean']] <- 0
    outputlist[['MSE']][['best']] <- 0
    outputlist[['genes']][['initial']] <- init.fit[['good.genes']]
    outputlist[['genes']][['mean']] <- init.fit[['good.genes']]
    outputlist[['genes']][['best']] <- init.fit[['good.genes']]
    
    return(outputlist)
  }
  
  
  outputlist[['genes']][['init']] <- init.fit[["good.genes"]]
  outputlist[['mse']][['init']] <- init.fit[["mse"]]
  
  
  
  best.model <- c()
  best.model$mse <- init.fit[["mse"]]
  outputlist[["MSE"]][['initial']] <- init.fit[["mse"]]
  
  # Finding the combination of genes which minimizes the MSE
  thresholds <- seq(35,60)
  t_best <- NA
  
  for (t in thresholds){
    i <- 0 #rounds of iterations
    #resets for new threshold
    genes <- names(MSPlist[[id]][1:n.genes, 1])
    MSE.old <- best.model[['mse']]
    MSE <- ""
    new_genes_assigned <- F
    
    #carries out the initial fit
    init.genes <- rank.mat(MSPlist, id = id,
                           gene.names = genes,
                           step = "mean",
                           phase = "initial",
                           threshold = t,
                           t.start = 0,
                           t.end = n.genes,
                           mapped = mapped)
    
    # iteratively replace genes until the MSE does not decrease (with an interval for MSE allowed)
    while ((MSE < (MSE.old * 1.01)) & (! MSE %in% mse[1:length(mse) - 1])){ # if the mse is seen before, break
      #actually does the fit, test and stuff
      
      #If there has been new genes assigned, then you can re-perform an initial fit
      if(new_genes_assigned==T){
        init.genes <- rank.mat(MSPlist, id,
                               gene.names = genes,
                               "mean",
                               "initial",
                               t,
                               0,
                               n.genes,mapped = mapped)
      }
      
      if (length(init.genes$good.genes) < 10){ # if there is below 10 good genes it is impossible to fit model and there is no reason to keep going, hence the break
        MSE <- init.genes$mse
        break
      }
      
      #if all the genes meet the threshold set up before then we are done
      #if we have 100 good genes then there is no reason to keep going
      if (length(init.genes$good.genes) == n.genes){
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- init.genes$good.genes
        
        break
      }
      
      # if this is the first iteration, then set the MSE.old as the current MSE
      if (i == 0) {
        MSE.old <- init.genes$mse
        mse.before[id] <- init.genes$mse
        
        # if not, then update the old MSE as the MSE of the from the previous iteration
      } else {
        MSE.old <- rotation$mse
      }
      
      rotation <- rank.mat(MSPlist, id = id,
                           gene.names = init.genes$good.genes,
                           step = "mean",
                           phase = "rotation",
                           threshold = t,
                           t.start = n.genes,
                           t.end = min(500, length(MSPlist[[id]][, 1])),
                           n.replace = (n.genes - length(init.genes$good.genes)),
                           mapped = mapped)
      
      # if the initial gene set is better than the rotation, then why continue? <You started out with the best thing
      if (is.na(rotation$mse) | (init.genes$mse < rotation$mse & i == 0)){
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- genes
        break
      }
      MSE <- rotation$mse
      rotation[['threshold']] <- t
      
      i <- i + 1
      mse <- c(mse, MSE.old)
      
      # Keeping the best model across iterations and thresholds if the new model is better than the best model across thresholds
      if ((MSE <= MSE.old) & (MSE < best.model$mse)){
        tmp.best.genes <- c(init.genes$good.genes, names(rotation$good.genes))
        if (length(tmp.best.genes)<100){
          print("There is not enough genes to rotate")
          break
        }
        best.model <- rotation
        best.genes <- tmp.best.genes
        genes <- best.genes
        t_best <- t
        #if you assign new genes, you can do another initial fit
        new_genes_assigned <- T
      }else{
        new_genes_assigned <- F
      }
    }#while over
  }
  
  outputlist[["genes"]][['mean']] <- best.genes
  outputlist[['MSE']][['mean']] <- best.model[['mse']]
  
  best.genes.step.one <- best.genes
  
  
  ############################################################## The 95-percentile gene refinement ##########################################################3
  # In the case that some genes perform extremely bad in a few samples, they will still affect the performance, and we want them replaced.
  
  best.old <- best.genes
  best.old.mse <- best.model$mse
  thres <- c(90,91,92,93,94,95,96,97,98) # The thresholds for 95-percentile ranks
  
  #if you never found a new better set, then you can just start over
  if(is.null(best.old)){
    best.old <- rank.mat(MSPlist, id,
                         gene.names = genes,
                         "mean",
                         "initial",
                         999,
                         0,
                         n.genes,mapped = mapped)[['good.genes']]
  }
  
  t_best <- NA
  for (t in thres){
    
    MSE <- 0
    MSE.old <- best.old.mse
    k <- 0
    new.genes <- best.old
    
    #initially, evaulated the init.quant
    genes_revised <- T 
    
    #while the MSE is improving (and has not been seen before)
    while((MSE < (MSE.old * 1.05)) & (! MSE %in% mse[1:length(mse) - 1])){
      
      if(genes_revised==T){
        init.quant <- rank.mat(MSPlist, id,new.genes, "quantile", "initial", t, 0, n.genes,mapped=mapped)
      }
      
      if (k != 0){
        MSE.old <- rot.quant$mse
      }else {
        MSE.old <- init.quant$mse
      }
      #redoes the fit with new genes found during rotation, if there had been found new during rotation
      if(genes_revised==T){
        init.quant <- rank.mat(MSPlist, id,new.genes, "quantile", "initial", t, 0, n.genes,mapped=mapped)
      }
      
      #if the initial genes are all included by the threshold, then we are good and break
      if (length(init.quant$good.genes) == n.genes){
        
        MSE <- init.quant$mse
        if (MSE < best.model$mse){
          
          best.model <- init.quant
          best.genes <- init.quant$good.genes}
        break
      }
      
      # if there is only 5 good genes it is impossible to fit model
      if (length(init.quant$good.genes) < 5){ # if there is only 5 good gene it is impossible to fit model
        if (init.quant$mse < best.model$mse){
          best.model <- init.quant
          best.genes <- new.genes
        }
        break
      }
      
      #rotates the genes
      rot.quant <- rank.mat(MSPlist, id,
                            init.quant$good.genes,
                            "quantile",
                            "rotation",
                            t,
                            n.genes,
                            min(500, length(MSPlist[[id]][, 1])),
                            (n.genes - length(init.quant$good.genes)),
                            mapped = mapped)
      new.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
      MSE <- rot.quant$mse
      k <- k + 1
      
      if (is.na(MSE) | (MSE > best.model$mse)){
        genes_revised <- F
        MSE <- best.model$mse
        mse <- c(mse, MSE)
      }else{
        mse <- c(mse, MSE)
        best.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
        if (length(best.genes)<100){
          print("There is not enough genes to rotate")
          break
        }
        best.model <- rot.quant
        genes_revised <- T
        t_best <- t
        
      }
    }
  }
  
  # Storing the best model
  outputlist[["genes"]][['best']] <- best.genes
  outputlist[["MSE"]][['best']] <- best.model$mse
  
  
  i <- match(id,ids)
  print(round(i/length(ids)*100,2))
  return(outputlist)
}
