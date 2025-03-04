# Initialize all per-cell-type comparisons.
rocAucDirSCORPION <- NULL
rocAucDirSCENIC <- NULL
ranges_all_file <- NULL
jaccard_plot <- NULL
indegree_plot <- NULL
outdegree_plot <- NULL

# Macrophages
sourceFileListMac <- c("1188_Ru1170gFreeze-sort_IGO_10034_8_dense.csv", 
                    rep("1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv", 2),
                    rep("1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv", 3), 
                    "RU661_TUMOR_dense.csv",
                    rep("RU682_NORMAL_dense.csv", 2), 
                    rep("RU682_TUMOR_dense.csv", 3),
                    "RU684_NORMAL_dense.csv", "Ru1135_dense.csv", "Ru1137_dense.csv")
ingroupNameMac <- rep("Macrophages", 15)
outgroupNameMac <- c(rep("Endothelial_cells", 2), "NK_cells", "B_cells", "Fibroblasts",
                     "T_cells", "T_cells", "NK_cells", "T_cells", "Endothelial_cells",
                     "Fibroblasts", "T_cells", "T_cells", "Fibroblasts", "Endothelial_cells")

# Endothelial cells
sourceFileListEndo <- c("1188_Ru1170gFreeze-sort_IGO_10034_8_dense.csv", 
                       rep("1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv", 2),
                       rep("RU682_TUMOR_dense.csv", 3),
                       "Ru1137_dense.csv")
ingroupNameEndo <- rep("Endothelial_cells", 7)
outgroupNameEndo <- c("Macrophages", "Macrophages", "NK_cells", "Fibroblasts",
                      "Macrophages", "T_cells", "Macrophages")

# NK cells
sourceFileListNK <- c(rep("1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv", 2),
                       rep("RU682_NORMAL_dense.csv", 2))
ingroupNameNK <- rep("NK_cells", 4)
outgroupNameNK <- c("Macrophages", "Endothelial_cells", "Macrophages", "T_cells")

# Fibroblasts
sourceFileListFib <- c(rep("1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv", 3), 
                       rep("RU682_TUMOR_dense.csv", 3),
                       "Ru1135_dense.csv")
ingroupNameFib <- rep("Fibroblasts", 7)
outgroupNameFib <- c("B_cells", "Macrophages", "T_cells", "Endothelial_cells",
                     "Macrophages", "T_cells", "Macrophages")

# T cells
sourceFileListT <- c(rep("1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv", 3), 
                       "RU661_TUMOR_dense.csv",
                       rep("RU682_NORMAL_dense.csv", 2), 
                       rep("RU682_TUMOR_dense.csv", 3),
                       "RU684_NORMAL_dense.csv")
ingroupNameT <- rep("T_cells", 10)
outgroupNameT <- c("B_cells", "Fibroblasts", "Macrophages", "Macrophages", "Macrophages",
                   "NK_cells", "Endothelial_cells", "Fibroblasts", "Macrophages",
                   "Macrophages")

# Read all AUC/ROC results from SCORPION
rocAucDirSCORPION <- "/Users/tae771/Documents/postdoc/single_cell_HTAN_MSK/SCORPION_roc_auc/"
# Macrophages
ingroupAUC_SCORPION_Mac <- lapply(1:length(ingroupNameMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListMac[i], "_", ingroupNameMac[i], 
                         "_", outgroupNameMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_macrophage.pdf"))
consolidatedAUCROC_SCORPION_Mac <- ConsolidateRobustness(ingroupAUC_SCORPION_Mac, xlab = "Macrophages to Others", ylab = "Macrophages to Macrophages",
                                                      minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_Mac <- GetResultRanges(ingroupAUC_SCORPION_Mac)
# Endothelial cells
ingroupAUC_SCORPION_Endo <- lapply(1:length(ingroupNameEndo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListEndo[i], "_", ingroupNameEndo[i], 
                                  "_", outgroupNameEndo[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_Endothelial_cells.pdf"))
consolidatedAUCROC_SCORPION_Endo <- ConsolidateRobustness(ingroupAUC_SCORPION_Endo, xlab = "Endothelial Cells to Others", ylab = "Endothelial Cells to Endothelial Cells",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_Endo <- GetResultRanges(ingroupAUC_SCORPION_Endo)
# NK cells
ingroupAUC_SCORPION_NK <- lapply(1:length(ingroupNameNK), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListNK[i], "_", ingroupNameNK[i], 
                                  "_", outgroupNameNK[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_NK_cells.pdf"))
consolidatedAUCROC_SCORPION_NK <- ConsolidateRobustness(ingroupAUC_SCORPION_NK, xlab = "NK Cells to Others", ylab = "NK Cells to NK Cells",
                                                        minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_NK <- GetResultRanges(ingroupAUC_SCORPION_NK)
# Fibroblasts
ingroupAUC_SCORPION_Fib <- lapply(1:length(ingroupNameFib), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListFib[i], "_", ingroupNameFib[i], 
                                  "_", outgroupNameFib[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_Fibroblasts.pdf"))
consolidatedAUCROC_SCORPION_Fib <- ConsolidateRobustness(ingroupAUC_SCORPION_Fib, xlab = "Fibroblasts to Others", ylab = "Fibroblasts to Fibroblasts",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_Fib <- GetResultRanges(ingroupAUC_SCORPION_Fib)
# T cells
ingroupAUC_SCORPION_T <- lapply(1:length(ingroupNameT), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListT[i], "_", ingroupNameT[i], 
                                  "_", outgroupNameT[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_T_cells.pdf"))
consolidatedAUCROC_SCORPION_T <- ConsolidateRobustness(ingroupAUC_SCORPION_T, xlab = "T Cells to Others", ylab = "T Cells to T Cells",
                                                       minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_T <- GetResultRanges(ingroupAUC_SCORPION_T)
# B cells
ingroupAUC_SCORPION_B <- lapply(1:length(ingroupNameT), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCORPION, sourceFileListT[i], "_", ingroupNameT[i], 
                                  "_", outgroupNameT[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCORPION, "consolidated_B_cells.pdf"))
consolidatedAUCROC_SCORPION_B <- ConsolidateRobustness(ingroupAUC_SCORPION_B, xlab = "B Cells to Others", ylab = "B Cells to B Cells",
                                                       minTextDistAsPercentage = 20)
dev.off()
ranges_SCORPION_B <- GetResultRanges(ingroupAUC_SCORPION_B)
# Consolidate
ranges_SCORPION <- do.call(rbind,  list(ranges_SCORPION_Mac, ranges_SCORPION_Endo, ranges_SCORPION_Fib,
                                        ranges_SCORPION_NK, ranges_SCORPION_T, ranges_SCORPION_B))
ranges_SCORPION$CellType <- c(rep("Macrophage", 6), rep("Endothelial Cells", 6),
                              rep("Fibroblasts", 6), rep("NK Cells", 6),
                              rep("T Cells", 6), rep("B Cells", 6))
ranges_SCORPION$Method <- rep("SCORPION", 36)
ranges_SCORPION <- ranges_SCORPION[,c("Method", "CellType", colnames(ranges_SCORPION)[1:(ncol(ranges_SCORPION)-2)])]

# Read all AUC/ROC results from SCENIC
rocAucDirSCENIC <- "/Users/tae771/Documents/postdoc/single_cell_HTAN_MSK/SCENIC_roc_auc/"
# Macrophages
ingroupAUC_SCENIC_Mac <- lapply(1:length(ingroupNameMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListMac[i], "_", ingroupNameMac[i], 
                                  "_", outgroupNameMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_macrophage.pdf"))
consolidatedAUCROC_SCENIC_Mac <- ConsolidateRobustness(ingroupAUC_SCENIC_Mac, xlab = "Macrophages to Others", ylab = "Macrophages to Macrophages",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_Mac <- GetResultRanges(ingroupAUC_SCENIC_Mac)
# Endothelial cells
ingroupAUC_SCENIC_Endo <- lapply(1:length(ingroupNameEndo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListEndo[i], "_", ingroupNameEndo[i], 
                                  "_", outgroupNameEndo[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_Endothelial_cells.pdf"))
consolidatedAUCROC_SCENIC_Endo <- ConsolidateRobustness(ingroupAUC_SCENIC_Endo, xlab = "Endothelial Cells to Others", ylab = "Endothelial Cells to Endothelial Cells",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_Endo <- GetResultRanges(ingroupAUC_SCENIC_Endo)
# NK cells
ingroupAUC_SCENIC_NK <- lapply(1:length(ingroupNameNK), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListNK[i], "_", ingroupNameNK[i], 
                                  "_", outgroupNameNK[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_NK_cells.pdf"))
consolidatedAUCROC_SCENIC_NK <- ConsolidateRobustness(ingroupAUC_SCENIC_NK, xlab = "NK Cells to Others", ylab = "NK Cells to NK Cells",
                                                        minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_NK <- GetResultRanges(ingroupAUC_SCENIC_NK)
# Fibroblasts
ingroupAUC_SCENIC_Fib <- lapply(1:length(ingroupNameFib), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListFib[i], "_", ingroupNameFib[i], 
                                  "_", outgroupNameFib[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_Fibroblasts.pdf"))
consolidatedAUCROC_SCENIC_Fib <- ConsolidateRobustness(ingroupAUC_SCENIC_Fib, xlab = "Fibroblasts to Others", ylab = "Fibroblasts to Fibroblasts",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_Fib <- GetResultRanges(ingroupAUC_SCENIC_Fib)
# T cells
ingroupAUC_SCENIC_T <- lapply(1:length(ingroupNameT), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListT[i], "_", ingroupNameT[i], 
                                  "_", outgroupNameT[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_T_cells.pdf"))
consolidatedAUCROC_SCENIC_T <- ConsolidateRobustness(ingroupAUC_SCENIC_T, xlab = "T Cells to Others", ylab = "T Cells to T Cells",
                                                       minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_T <- GetResultRanges(ingroupAUC_SCENIC_T)
# B cells
ingroupAUC_SCENIC_B <- lapply(1:length(ingroupNameT), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSCENIC, sourceFileListT[i], "_", ingroupNameT[i], 
                                  "_", outgroupNameT[i], "_results.csv")))
})
pdf(paste0(rocAucDirSCENIC, "consolidated_B_cells.pdf"))
consolidatedAUCROC_SCENIC_B <- ConsolidateRobustness(ingroupAUC_SCENIC_B, xlab = "B Cells to Others", ylab = "B Cells to B Cells",
                                                     minTextDistAsPercentage = 20)
dev.off()
ranges_SCENIC_B <- GetResultRanges(ingroupAUC_SCENIC_B)
# Consolidate
ranges_SCENIC <- do.call(rbind,  list(ranges_SCENIC_Mac, ranges_SCENIC_Endo, ranges_SCENIC_Fib,
                                      ranges_SCENIC_NK, ranges_SCENIC_T, ranges_SCENIC_B))
ranges_SCENIC$CellType <- c(rep("Macrophage", 6), rep("Endothelial Cells", 6),
                            rep("Fibroblasts", 6), rep("NK Cells", 6),
                            rep("T Cells", 6), rep("B Cells", 6))
ranges_SCENIC$Method <- rep("SCENIC", 36)
ranges_SCENIC <- ranges_SCENIC[,c("Method", "CellType", colnames(ranges_SCENIC)[1:(ncol(ranges_SCENIC)-2)])]

# Consolidate all files.
ranges_all <- do.call(rbind, list(ranges_SCORPION, ranges_SCENIC))
write.csv(ranges_all, ranges_all_file)
ranges_all <- read.csv(ranges_all_file)

# Write a function to make bar plots from consolidated files.
pdf(jaccard_plot)
MakePerformanceBarPlot(ranges_all, "Jaccard")
dev.off()
pdf(indegree_plot)
MakePerformanceBarPlot(ranges_all, "InDegree")
dev.off()
pdf(outdegree_plot)
MakePerformanceBarPlot(ranges_all, "OutDegree")
dev.off()