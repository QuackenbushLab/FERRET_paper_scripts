library("FERRET")

# Variables
scGeneRai_results <- NULL
scGeneRai_results_formatted <- NULL
GRISLI_results <- NULL
scSGL_results <- NULL
PIDC_results <- NULL
LEAP_results <- NULL
SINGE_results <- NULL
roc_auc_parent_dir <- NULL

# Add source and target columns.
dir.create(scGeneRai_results_formatted)
for(file in list.files(scGeneRai_results)){
  dat <- utils::read.csv(paste0(scGeneRai_results, file),
                  row.names = 1)
  sources <- unlist(lapply(rownames(dat), function(row){
    return(strsplit(row, "_")[[1]][1])
  }))
  targets <- unlist(lapply(rownames(dat), function(row){
    return(strsplit(row, "_")[[1]][2])
  }))
  newDf <- data.frame(source = sources, target = targets, score = dat$LRP)
  str(newDf)
  rownames(newDf) <- rownames(dat)
  write.csv(newDf, paste0(scGeneRai_results_formatted, file))
}

# Convert matrix into GRN. Here, we extract the matrix at L = 20.
dir.create(GRILSI_results)
for(file in list.files(GRISLI_results)){
  dat <- utils::read.csv(paste0(GRISLI_results, file), header = FALSE)[,1:100]
  colnames(dat) <- paste0("PC", 1:100)
  rownames(dat) <- paste0("PC", 1:100)
  dat2 <- reshape2::melt(dat)
  dat2$regulator <- as.character(dat2$variable)
  dat2$target <- unlist(lapply(1:100, function(i){return(paste0("PC", 1:100))}))
  dat2 <- dat2[,c("regulator", "target", "value")]
  str(dat2)
  write.csv(dat2, paste0(GRILSI_results, file, ".csv"), row.names = FALSE)
}

# Load data.
networksScSGL <- LoadResults(resultDirectory = scSGL_results, 
                        interpretationOfNegative = "inhibitory")
networksPIDC <- LoadResults(resultDirectory = PIDC_results, 
                        interpretationOfNegative = "inhibitory")
networksLEAP <- LoadResults(resultDirectory = LEAP_results, 
                        interpretationOfNegative = "inhibitory")
networksScGeneRai <- LoadResults(resultDirectory = scGeneRai_results_formatted, 
                        interpretationOfNegative = "inhibitory")
networksSINGE <- LoadResults(resultDirectory = SINGE_results, 
                        interpretationOfNegative = "inhibitory", isTabDelimited = TRUE, firstColumnIsRowname = FALSE)
networksGRISLI <- LoadResults(resultDirectory = GRILSI_results, 
                        interpretationOfNegative = "inhibitory", isTabDelimited = FALSE, firstColumnIsRowname = FALSE)

# Create a repository for storing results.
roc_auc_dir_scSGL <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/scSGL_results_roc_auc_updatedFERRET/"
roc_auc_dir_PIDC <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/PIDC_results_roc_auc_updatedFERRET/"
roc_auc_dir_LEAP <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/LEAP_results_roc_auc_updatedFERRET/"
roc_auc_dir_scGeneRai <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/scGeneRai_results_roc_auc_updatedFERRET/"
roc_auc_dir_SINGE <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/SINGE_results_roc_auc_updatedFERRET/"
roc_auc_dir_GRISLI <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc_updatedFERRET/"

if(!file.exists(roc_auc_dir)){
  dir.create(roc_auc_dir)
}

# Helper function to build comparison objects.
getFileNames <- function(rootName, iterationCount, foldCount, ending, networks){
  comparisonList <- unlist(lapply(1:iterationCount, function(i){
    return(lapply(1:foldCount, function(j){
      retVal <- paste(rootName, i, j, ending, sep = "_")
      if(!file.exists(paste0(networks@directory, "/", retVal))){
        retVal <- NULL
      }
      return(retVal)
    }))
  }))
  return(comparisonList)
}

# Run all comparisons.
sourceFileList <- c(rep("C3L-00359-01", 2),
                    rep("C3N-00662-03", 2),
                    rep("C3N-01814-01", 10),
                    rep("C3N-01904-02", 2),
                    rep("C3N-02181-02", 10),
                    rep("C3N-02784-01", 2))
ingroupName <- c("Macrophages", "Podocytes",
                 "Oligodendrocytes", "Neurons",
                 rep("Interneurons", 2),
                 rep("Macrophages", 3),
                 rep("Neurons", 2),
                 rep("Oligodendrocytes", 3),
                 "Macrophages", "Podocytes",
                 rep("Interneurons", 2),
                 rep("Macrophages", 3),
                 rep("Neurons", 2),
                 rep("Oligodendrocytes", 3),
                 "Macrophages", "Oligodendrocytes")
outgroupName <- c("Podocytes", "Macrophages",
                  "Neurons", "Oligodendrocytes",
                  "Macrophages", "Oligodendrocytes",
                  "Interneurons", "Neurons", "Oligodendrocytes",
                  "Macrophages", "Oligodendrocytes",
                  "Macrophages", "Neurons", "Interneurons",
                  "Podocytes", "Macrophages",
                  "Macrophages", "Oligodendrocytes",
                  "Interneurons", "Neurons", "Oligodendrocytes",
                  "Macrophages", "Oligodendrocytes",
                  "Macrophages", "Neurons", "Interneurons",
                  "Oligodendrocytes", "Macrophages")
uniqueSamples <- unique(sourceFileList)

# Do AUC.
for(i in 1:length(sourceFileList)){
    tryCatch({
       #Compute the AUC and plot the curve.
      comparisonsScSGL <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv.csv", sep = "_"),
                                          ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv.csv", networksScSGL),
                                          outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv.csv", networksScSGL),
                                          results = networks)
      comparisonsPIDC <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv", sep = "_"),
                                            ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv", networksPIDC),
                                            outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv", networksPIDC),
                                           results = networks)
      comparisonsLEAP <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv", sep = "_"),
                                           ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv", networksLEAP),
                                           outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv", networksLEAP),
                                           results = networks)
      comparisonsScGeneRai <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv_consolidated.csv", sep = "_"),
                                             ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv_consolidated.csv", networksScGeneRai),
                                             outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv_consolidated.csv", networksScGeneRai),
                                             results = networks)
      comparisonsSINGE <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv.txt", sep = "_"),
                                          ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv.txt", networksSINGE),
                                          outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv.txt", networksSINGE),
                                          results = networks)
      comparisonsGRISLI <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], ingroupName[i], "original_pca.tsv.csv", sep = "_"),
                                         ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], sep = "_"), 5, 10, "pca.tsv.csv", networksGRISLI),
                                         outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], sep = "_"), 5, 10, "pca.tsv.csv", networksGRISLI),
                                         results = networks)
      pdf(paste0(roc_auc_dir_scSGL, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksScSGL,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                 ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_scSGL, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_scSGL, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_PIDC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksPIDC,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                     ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_PIDC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_PIDC, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_LEAP, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksLEAP,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                     ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_LEAP, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_LEAP, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_scGeneRai, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksscGeneRai,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                     ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_scGeneRai, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_scGeneRai, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_SINGE, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksSINGE,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                     ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_SINGE, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_SINGE, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      pdf(paste0(roc_auc_dir_GRISLI, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), ".pdf")))
      result <- ComputeRobustnessAUC(results = networksGRISLI,comparisons = comparisons,xlab = paste(ingroupName[i], "to", outgroupName[i]),
                                     ylab = paste(ingroupName[i], "to", ingroupName[i]), metric = c("jaccard", "in-degree", "out-degree"))
      dev.off()
      WriteRobustnessAUC(result, paste0(roc_auc_dir_GRISLI, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      print(paste0(roc_auc_dir_GRISLI, paste0("/", paste(sourceFileList[i], ingroupName[i], outgroupName[i], sep = "_"), "_results.csv")))
      
    }, error = function(cond){
      print(paste("Skipping iteration", i))
    })
}