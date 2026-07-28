outputDir <- NULL

# Read all AUC/ROC results for random networks.
allResults <- do.call(c, lapply(1:5, function(i){
  return(lapply(1:20, function(j){
    return(ReadRobustnessAUC(paste0(outputDir, "/roc_", i, "_", j, "_results.csv")))
  }))
}))
pdf(paste0(outputDir, "/results_consolidated.pdf"))
ConsolidateRobustness(allResults, xlab = "Out-Group Similarity", ylab = "In-Group Similarity",
                                                      minTextDistAsPercentage = 20)
dev.off()