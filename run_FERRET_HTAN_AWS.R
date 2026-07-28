library("FERRET")

# Load data.
networks_SCORPION_file <- NULL
networks_SCENIC_file <- NULL
networksSCORPION <- LoadResults(resultDirectory = networks_SCORPION_file, 
                        interpretationOfNegative = "poor", firstColumnIsRowname = TRUE)
networksSCENIC <- LoadResults(resultDirectory = networks_SCENIC_file, 
                        interpretationOfNegative = "inhibitory", firstColumnIsRowname = TRUE)

# Create a repository for storing results.
roc_auc_dir_SCORPION <- NULL
roc_auc_dir_SCENIC <- NULL

if(!file.exists(roc_auc_dir)){
  dir.create(roc_auc_dir)
}

# Helper function to build comparison objects.
foldCount <- 5
getFileNames <- function(rootName, iterationCount, foldCount, ending){
  comparisonList <- unlist(lapply(1:foldCount, function(j){
    retVal <- paste(rootName, 20, j, ending, sep = "_")
    if(!file.exists(paste0(networks@directory, "/", retVal))){
      retVal <- NULL
    }
    return(retVal)
  }))
  return(comparisonList)
}

# Run all comparisons.
sourceFileList <- c(rep("1188_Ru1170gFreeze-sort_IGO_10034_8_dense.csv", 2),
                    rep("1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv", 6),
                    rep("1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv", 12),
                    rep("RU661_TUMOR_dense.csv", 2),
                    rep("RU682_NORMAL_dense.csv", 6),
                    rep("RU682_TUMOR_dense.csv", 12),
                    rep("RU684_NORMAL_dense.csv", 2),
                    rep("Ru1135_dense.csv", 2),
                    rep("Ru1137_dense.csv", 2))
ingroupName <- c("Endothelial_cells", "Macrophages",
                 rep("Endothelial_cells", 2),
                 rep("Macrophages", 2),
                 rep("NK_cells", 2),
                 rep("B_cells", 3),
                 rep("Fibroblasts", 3),
                 rep("Macrophages", 3),
                 rep("T_cells", 3),
                 "Macrophages", "T_cells",
                 rep("Macrophages", 2),
                 rep("NK_cells", 2),
                 rep("T_cells", 2),
                 rep("Endothelial_cells", 3),
                 rep("Fibroblasts", 3),
                 rep("Macrophages", 3),
                 rep("T_cells", 3),
                 "Macrophages", "T_cells",
                 "Fibroblasts", "Macrophages",
                 "Endothelial_cells", "Macrophages")
outgroupName <- c("Macrophages", "Endothelial_cells",
                  "Macrophages", "NK_cells",
                  "Endothelial_cells", "NK_cells",
                  "Endothelial_cells", "Macrophages",
                  "Fibroblasts", "Macrophages", "T_cells",
                  "B_cells", "Macrophages", "T_cells",
                  "B_cells", "Fibroblasts", "T_cells",
                  "B_cells", "Fibroblasts", "Macrophages",
                  "T_cells", "Macrophages",
                  "NK_cells", "T_cells",
                  "Macrophages", "T_cells",
                  "Macrophages", "NK_cells",
                  "Fibroblasts", "Macrophages", "T_cells",
                  "Endothelial_cells", "Macrophages", "T_cells",
                  "Endothelial_cells", "Fibroblasts", "T_cells",
                  "Endothelial_cells", "Fibroblasts", "Macrophages",
                  "T_cells", "Macrophages",
                  "Macrophages", "Fibroblasts",
                  "Macrophages", "Endothelial_cells")
uniqueSamples <- unique(sourceFileList)

# Do AUC.
for(i in 1:length(sourceFileList)){
      # Compute the AUC and plot the curve.
      comparisons <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "logcounts.csv_original.tsv.csv", sep = "_"),
                                          ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv.csv"),
                                          outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv.csv"),
                                          results = networks)
      pdf(paste0(roc_auc_dir_SCORPION, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      resultSCORPION <- ComputeRobustnessAUC(results = networksSCORPION,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                 ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_SCORPION, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_SCORPION, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_SCENIC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      resultSCENIC <- ComputeRobustnessAUC(results = networksSCENIC,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                           ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_SCENIC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_SCENIC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
}