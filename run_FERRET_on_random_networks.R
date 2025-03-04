# First, generate 100 random scale-free networks using the Barabasi-Albert method
inputDir <- NULL
outputDir <- NULL
if(!file.exists(inputDir)){
  dir.create(inputDir)
}
if(!file.exists(outputDir)){
  dir.create(outputDir)
}
if(!require(igraph)){
  install.packages("igraph")
}
library("igraph")

# Generate the networks.
set.seed(1)
for(i in 1:100){
  
  # Build network.
  network <- sample_pa(n = 1000, directed = TRUE, algorithm = "psumtree")

  # Format as data frame.
  df <- as_long_data_frame(network)
  colnames(df) <- c("tf", "gene")
  
  # Add scores.
  df$score <- runif(n = nrow(df), min = 0.1, max = 1)
  
  # Save data frame.
  write.csv(df, paste0(inputDir, "/", i, ".csv"))
}

# Load the networks for FERRET.
library("FERRET")
networks <- LoadResults(resultDirectory = inputDir, interpretationOfNegative = "inhibitory")

# Separate into folds.
shuffledNets <- networks@results[sample(length(networks@results))]
totalFoldCount <- 5
folds <- lapply(1:totalFoldCount, function(i){
  step <- length(shuffledNets) / totalFoldCount
  start <- (i - 1) * step + 1
  end <- i * step
  print(paste(start, end))
  return(shuffledNets[start:end])
})

# For each fold, run FERRET including each network as the source.
for(j in 1:length(folds)){
  fold <- folds[[j]]
  
  # Loop through each network in the fold.
  for(i in 1:length(fold)){
    
    # Include all networks not in the fold in the ingroup and outgroup.
    network <- names(fold)[i]
    everythingElse <- setdiff(names(fold), network)
    
    # Create the comparisons.
    comparisons <- BuildComparisonObject(sourceNetwork = network, ingroupToCompare = everythingElse[1:9],
                                         outgroupToCompare = everythingElse[10:19], results = networks)
    
    # Run FERRET.
    pdf(paste0(outputDir, "/roc_", j, "_", i, ".pdf"))
    result <- ComputeRobustnessAUC(results = networks,comparisons = comparisons,xlab = "Ingroup to Outgroup",
                                   ylab = "Ingroup to Ingroup", metric = c("jaccard", "in-degree", "out-degree"))
    dev.off()
    WriteRobustnessAUC(result, paste0(outputDir, "/roc_", j, "_", i, "_results.csv"))
    
  }
}