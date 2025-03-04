# Initialize all comparisons.
sourceFileListRCC <- c(rep("C3L-00359-01", 2), 
                    rep("C3N-01904-02", 2))
ingroupNameRCC <- c("Macrophages", "Podocytes", 
                 "Macrophages", "Podocytes")
outgroupNameRCC <- c("Podocytes", "Macrophages", 
                  "Podocytes", "Macrophages")
sourceFileListGBM <- c(rep("C3N-00662-03", 2), 
                       rep("C3N-01814-01", 10),
                       rep("C3N-02181-02", 10),
                       rep("C3N-02784-01", 2))
ingroupNameGBM <- c("Oligodendrocytes", "Neurons", 
                    rep("Interneurons", 2),
                    rep("Macrophages", 3),
                    rep("Neurons", 2),
                    rep("Oligodendrocytes", 3),
                    rep("Interneurons", 2),
                    rep("Macrophages", 3),
                    rep("Neurons", 2),
                    rep("Oligodendrocytes", 3),
                    "Macrophages", "Oligodendrocytes")
outgroupNameGBM <- c("Neurons", "Oligodendrocytes", 
                     "Macrophages", "Oligodendrocytes",
                     "Interneurons", "Neurons", "Oligodendrocytes",
                     "Macrophages", "Oligodendrocytes",
                     "Macrophages", "Neurons", "Interneurons",
                     "Macrophages", "Oligodendrocytes",
                     "Interneurons", "Neurons", "Oligodendrocytes",
                     "Macrophages", "Oligodendrocytes",
                     "Macrophages", "Neurons", "Interneurons",
                     "Oligodendrocytes", "Macrophages")


# Read all AUC/ROC results from scSGL.
rocAucDirScSGL <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/scSGL_results_roc_auc_updatedFERRET/"
ingroupAUC_scSGL_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                         "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/scSGL_results_roc_auc_subset/RCC_consolidated.pdf")
consolidatedAUCROC_scSGL_RCC <- ConsolidateRobustness(ingroupAUC_scSGL_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                      minTextDistAsPercentage = 20)
scSGL_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirScSGL, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "scSGL"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
dev.off()
ingroupAUC_scSGL_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/scSGL_results_roc_auc_subset/GBM_consolidated.pdf")
consolidatedAUCROC_scSGL_GBM <- ConsolidateRobustness(ingroupAUC_scSGL_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                      minTextDistAsPercentage = 20)
scSGL_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirScSGL, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "scSGL"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
dev.off()

# Read all AUC/ROC results from PIDC.
rocAucDirPIDC <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/PIDC_results_roc_auc_updatedFERRET/"
ingroupAUC_PIDC_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/PIDC_results_roc_auc/RCC_consolidated.pdf")
consolidatedAUCROC_PIDC_RCC <- ConsolidateRobustness(ingroupAUC_PIDC_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                     minTextDistAsPercentage = 20)
PIDC_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirPIDC, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                             "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "PIDC"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
dev.off()
ingroupAUC_PIDC_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/PIDC_results_roc_auc/GBM_consolidated.pdf")
consolidatedAUCROC_PIDC_GBM <- ConsolidateRobustness(ingroupAUC_PIDC_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                     minTextDistAsPercentage = 20)
dev.off()
PIDC_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirPIDC, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "PIDC"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))

# Read all AUC/ROC results from LEAP.
rocAucDirLEAP <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/LEAP_results_roc_auc_updatedFERRET/"
ingroupAUC_LEAP_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/LEAP_results_roc_auc/RCC_consolidated.pdf")
consolidatedAUCROC_LEAP_RCC <- ConsolidateRobustness(ingroupAUC_LEAP_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                     minTextDistAsPercentage = 20)
dev.off()
LEAP_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirLEAP, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                             "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "LEAP"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
ingroupAUC_LEAP_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/LEAP_results_roc_auc/GBM_consolidated.pdf")
consolidatedAUCROC_LEAP_GBM <- ConsolidateRobustness(ingroupAUC_LEAP_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                     minTextDistAsPercentage = 20)
dev.off()
LEAP_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirLEAP, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "LEAP"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))

# Read all AUC/ROC results from scGeneRai.
rocAucDirscGeneRai <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/scGeneRai_results_roc_auc_updatedFERRET/"
ingroupAUC_scGeneRai_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirscGeneRai, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/scGeneRai_results_roc_auc/RCC_consolidated.pdf")
consolidatedAUCROC_scGeneRai_RCC <- ConsolidateRobustness(ingroupAUC_scGeneRai_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
scGeneRai_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirscGeneRai, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                             "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "scGeneRai"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
ingroupAUC_scGeneRai_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirscGeneRai, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/scGeneRai_results_roc_auc/GBM_consolidated.pdf")
consolidatedAUCROC_scGeneRai_GBM <- ConsolidateRobustness(ingroupAUC_scGeneRai_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
scGeneRai_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirscGeneRai, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "scGeneRai"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))

# Read all AUC/ROC results from SINGE.
rocAucDirSINGE <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/SINGE_results_roc_auc_updatedFERRET/"
ingroupAUC_SINGE_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/SINGE_results_roc_auc/RCC_consolidated.pdf")
consolidatedAUCROC_scGeneRai_RCC <- ConsolidateRobustness(ingroupAUC_SINGE_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
SINGE_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirSINGE, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                             "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "SINGE"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
ingroupAUC_SINGE_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/SINGE_results_roc_auc/GBM_consolidated.pdf")
consolidatedAUCROC_SINGE_GBM <- ConsolidateRobustness(ingroupAUC_SINGE_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
SINGE_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirSINGE, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "SINGE"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))

# Read all AUC/ROC results from GRISLI
rocAucDirGRISLI <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc_updatedFERRET/"
ingroupAUC_GRISLI_RCC <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                                  "_", outgroupNameRCC[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc/RCC_consolidated.pdf")
consolidatedAUCROC_scGeneRai_RCC <- ConsolidateRobustness(ingroupAUC_GRISLI_RCC, xlab = "Out-Group Similarity (RCC)", ylab = "In-Group Similarity (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
GRISLI_RCC_df <- do.call(rbind, lapply(1:length(ingroupNameRCC), function(i){
  results <- read.csv(paste0(rocAucDirGRISLI, sourceFileListRCC[i], "_", ingroupNameRCC[i], 
                             "_", outgroupNameRCC[i], "_results.csv"))
  results$dataset <- "RCC"
  results$method <- "GRISLI"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))
ingroupAUC_GRISLI_GBM <- lapply(1:length(ingroupNameRCC), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                                  "_", outgroupNameGBM[i], "_results.csv")))
})
pdf("/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc/GBM_consolidated.pdf")
consolidatedAUCROC_GRISLI_GBM <- ConsolidateRobustness(ingroupAUC_GRISLI_GBM, xlab = "Out-Group Similarity (GBM)", ylab = "In-Group Similarity (GBM)",
                                                      minTextDistAsPercentage = 20)
dev.off()
GRISLI_GBM_df <- do.call(rbind, lapply(1:length(ingroupNameGBM), function(i){
  results <- read.csv(paste0(rocAucDirGRISLI, sourceFileListGBM[i], "_", ingroupNameGBM[i], 
                             "_", outgroupNameGBM[i], "_results.csv"))
  results$dataset <- "GBM"
  results$method <- "GRISLI"
  results <- results[,c("dataset", "method", "X", "Jaccard", "InDegree", "OutDegree", "cutoff", "inOrOut")]
}))

write.csv(do.call(rbind, list(scSGL_RCC_df, PIDC_RCC_df, LEAP_RCC_df, scGeneRai_RCC_df, SINGE_RCC_df, GRISLI_RCC_df)),
          "/Users/tae771/Documents/postdoc/single_cell_CPTAC/allResults_RCC.csv")
write.csv(do.call(rbind, list(scSGL_GBM_df, PIDC_GBM_df, LEAP_GBM_df, scGeneRai_GBM_df, SINGE_GBM_df, GRISLI_GBM_df)),
          "/Users/tae771/Documents/postdoc/single_cell_CPTAC/allResults_GBM.csv")
