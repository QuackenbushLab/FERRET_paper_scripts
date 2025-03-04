# Initialize all per-cell-type comparisons.
rocAucDirScSGL <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/scSGL_results_roc_auc_updatedFERRET/"
rocAucDirPIDC <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/PIDC_results_roc_auc_updatedFERRET/"
rocAucDirLEAP <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/LEAP_results_roc_auc_updatedFERRET/"
rocAucDirScGeneRai <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/ScGeneRai_results_roc_auc_updatedFERRET/"
rocAucDirSINGE <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/SINGE_results_roc_auc_updatedFERRET/"
rocAucDirGRISLI <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc_updatedFERRET/"
ranges_all_file <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/ranges_all.csv"
jaccard_plot <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/performanceBarPlot_jaccard.pdf"
indegree_plot <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/performanceBarPlot_jaccard.pdf"
outdegree_plot <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/performanceBarPlot_jaccard.pdf"

# RCC - macrophages and podocytes
sourceFileListRCC <- c("C3L-00359-01", "C3N-01904-02", 2)
ingroupNameRCCMac <- rep("Macrophages", 2)
outgroupNameRCCMac <- rep("Podocytes", 2)
ingroupNameRCCPod <- outgroupNameRCCMac
outgroupNameRCCPod <- ingroupNameRCCMac

# GBM, oligodendrocytes
sourceFileListGBMOligo <- c("C3N-00662-03", 
                       rep("C3N-01814-01", 3),
                       rep("C3N-02181-02", 3),
                       "C3N-02784-01")
ingroupNameGBMOligo <- c("Oligodendrocytes", 
                    rep("Oligodendrocytes", 3),
                    rep("Oligodendrocytes", 3))
outgroupNameGBMOligo <- c("Neurons",  
                     "Macrophages", "Neurons", "Interneurons",
                     "Macrophages", "Neurons", "Interneurons",
                     "Macrophages")

# GBM, neurons
sourceFileListGBMNeur <- c("C3N-00662-03", 
                       rep("C3N-01814-01", 2),
                       rep("C3N-02181-02", 2))
ingroupNameGBMNeur <- c("Neurons", 
                    rep("Neurons", 2),
                    rep("Neurons", 2))
outgroupNameGBMNeur <- c("Oligodendrocytes", 
                     "Macrophages", "Oligodendrocytes",
                     "Macrophages", "Oligodendrocytes")

# GBM, interneurons
sourceFileListGBMInt <- c(rep("C3N-01814-01", 2),
                       rep("C3N-02181-02", 2))
ingroupNameGBMInt <- c(rep("Interneurons", 2),
                    rep("Interneurons", 2))
outgroupNameGBMInt <- c("Macrophages", "Oligodendrocytes",
                     "Macrophages", "Oligodendrocytes")

# GBM, macrophages
sourceFileListGBMMac <- c(rep("C3N-01814-01", 3),
                       rep("C3N-02181-02", 3),
                       "C3N-02784-01")
ingroupNameGBMMac <- c(rep("Macrophages", 3),
                    rep("Macrophages", 3),
                    "Macrophages")
outgroupNameGBMMac <- c("Interneurons", "Neurons", "Oligodendrocytes",
                     "Interneurons", "Neurons", "Oligodendrocytes",
                     "Oligodendrocytes")


# Read all AUC/ROC results from scSGL.
# RCC Macrophages
ingroupAUC_scSGL_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                         "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_scSGL_RCC_Mac <- ConsolidateRobustness(ingroupAUC_scSGL_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                      minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_RCC_Mac <- GetResultRanges(ingroupAUC_scSGL_RCC_Mac)
# RCC Podocytes
ingroupAUC_scSGL_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_scSGL_RCC_Pod <- ConsolidateRobustness(ingroupAUC_scSGL_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_RCC_Pod <- GetResultRanges(ingroupAUC_scSGL_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_scSGL_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_scSGL_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_scSGL_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                            minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_GBM_Oligo <- GetResultRanges(ingroupAUC_scSGL_GBM_Oligo)
# GBM neurons
ingroupAUC_scSGL_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_scSGL_GBM_Neur <- ConsolidateRobustness(ingroupAUC_scSGL_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_GBM_Neur <- GetResultRanges(ingroupAUC_scSGL_GBM_Neur)
# GBM interneurons
ingroupAUC_scSGL_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_scSGL_GBM_Int <- ConsolidateRobustness(ingroupAUC_scSGL_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_GBM_Int <- GetResultRanges(ingroupAUC_scSGL_GBM_Int)
# GBM macrophages
ingroupAUC_scSGL_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScSGL, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirScSGL, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_scSGL_GBM_Mac <- ConsolidateRobustness(ingroupAUC_scSGL_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_scSGL_GBM_Mac <- GetResultRanges(ingroupAUC_scSGL_GBM_Mac)
# Consolidate results.
ranges_scSGL <- do.call(rbind,  list(ranges_scSGL_RCC_Mac, ranges_scSGL_RCC_Pod, ranges_scSGL_GBM_Oligo,
                              ranges_scSGL_GBM_Neur, ranges_scSGL_GBM_Int, ranges_scSGL_GBM_Mac))
ranges_scSGL$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                           rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                           rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_scSGL$Method <- rep("scSGL", 36)
ranges_scSGL <- ranges_scSGL[,c("Method", "CellType", colnames(ranges_scSGL)[1:(ncol(ranges_scSGL)-2)])]

# Read all AUC/ROC results from PIDC.
# RCC Macrophages
ingroupAUC_PIDC_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                                  "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_PIDC_RCC_Mac <- ConsolidateRobustness(ingroupAUC_PIDC_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_RCC_Mac <- GetResultRanges(ingroupAUC_PIDC_RCC_Mac)
# RCC Podocytes
ingroupAUC_PIDC_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_PIDC_RCC_Pod <- ConsolidateRobustness(ingroupAUC_PIDC_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_RCC_Pod <- GetResultRanges(ingroupAUC_PIDC_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_PIDC_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_PIDC_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_PIDC_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_GBM_Oligo <- GetResultRanges(ingroupAUC_PIDC_GBM_Oligo)
# GBM neurons
ingroupAUC_PIDC_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_PIDC_GBM_Neur <- ConsolidateRobustness(ingroupAUC_PIDC_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_GBM_Neur <- GetResultRanges(ingroupAUC_PIDC_GBM_Neur)
# GBM interneurons
ingroupAUC_PIDC_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_PIDC_GBM_Int <- ConsolidateRobustness(ingroupAUC_PIDC_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_GBM_Int <- GetResultRanges(ingroupAUC_PIDC_GBM_Int)
# GBM macrophages
ingroupAUC_PIDC_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirPIDC, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirPIDC, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_PIDC_GBM_Mac <- ConsolidateRobustness(ingroupAUC_PIDC_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_PIDC_GBM_Mac <- GetResultRanges(ingroupAUC_PIDC_GBM_Mac)
# Consolidate results.
ranges_PIDC <- do.call(rbind,  list(ranges_PIDC_RCC_Mac, ranges_PIDC_RCC_Pod, ranges_PIDC_GBM_Oligo,
                                    ranges_PIDC_GBM_Neur, ranges_PIDC_GBM_Int, ranges_PIDC_GBM_Mac))
ranges_PIDC$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                          rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                          rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_PIDC$Method <- rep("PIDC", 36)
ranges_PIDC <- ranges_PIDC[,c("Method", "CellType", colnames(ranges_PIDC)[1:(ncol(ranges_PIDC)-2)])]

# Read all AUC/ROC results from LEAP.
# RCC Macrophages
ingroupAUC_LEAP_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                                  "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_LEAP_RCC_Mac <- ConsolidateRobustness(ingroupAUC_LEAP_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_RCC_Mac <- GetResultRanges(ingroupAUC_LEAP_RCC_Mac)
# RCC Podocytes
ingroupAUC_LEAP_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_LEAP_RCC_Pod <- ConsolidateRobustness(ingroupAUC_LEAP_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_RCC_Pod <- GetResultRanges(ingroupAUC_LEAP_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_LEAP_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_LEAP_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_LEAP_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_GBM_Oligo <- GetResultRanges(ingroupAUC_LEAP_GBM_Oligo)
# GBM neurons
ingroupAUC_LEAP_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_LEAP_GBM_Neur <- ConsolidateRobustness(ingroupAUC_LEAP_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_GBM_Neur <- GetResultRanges(ingroupAUC_LEAP_GBM_Neur)
# GBM interneurons
ingroupAUC_LEAP_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_LEAP_GBM_Int <- ConsolidateRobustness(ingroupAUC_LEAP_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_GBM_Int <- GetResultRanges(ingroupAUC_LEAP_GBM_Int)
# GBM macrophages
ingroupAUC_LEAP_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirLEAP, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirLEAP, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_LEAP_GBM_Mac <- ConsolidateRobustness(ingroupAUC_LEAP_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                         minTextDistAsPercentage = 20)
dev.off()
ranges_LEAP_GBM_Mac <- GetResultRanges(ingroupAUC_LEAP_GBM_Mac)
# Consolidate results.
ranges_LEAP <- do.call(rbind,  list(ranges_LEAP_RCC_Mac, ranges_LEAP_RCC_Pod, ranges_LEAP_GBM_Oligo,
                                    ranges_LEAP_GBM_Neur, ranges_LEAP_GBM_Int, ranges_LEAP_GBM_Mac))
ranges_LEAP$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                          rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                          rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_LEAP$Method <- rep("LEAP", 36)
ranges_LEAP <- ranges_LEAP[,c("Method", "CellType", colnames(ranges_LEAP)[1:(ncol(ranges_LEAP)-2)])]


# Read all AUC/ROC results from ScGeneRai.
# RCC Macrophages
ingroupAUC_ScGeneRai_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                                  "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_ScGeneRai_RCC_Mac <- ConsolidateRobustness(ingroupAUC_ScGeneRai_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                              minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_RCC_Mac <- GetResultRanges(ingroupAUC_ScGeneRai_RCC_Mac)
# RCC Podocytes
ingroupAUC_ScGeneRai_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_ScGeneRai_RCC_Pod <- ConsolidateRobustness(ingroupAUC_ScGeneRai_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                              minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_RCC_Pod <- GetResultRanges(ingroupAUC_ScGeneRai_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_ScGeneRai_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_ScGeneRai_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_ScGeneRai_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                                minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_GBM_Oligo <- GetResultRanges(ingroupAUC_ScGeneRai_GBM_Oligo)
# GBM neurons
ingroupAUC_ScGeneRai_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_ScGeneRai_GBM_Neur <- ConsolidateRobustness(ingroupAUC_ScGeneRai_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                               minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_GBM_Neur <- GetResultRanges(ingroupAUC_ScGeneRai_GBM_Neur)
# GBM interneurons
ingroupAUC_ScGeneRai_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_ScGeneRai_GBM_Int <- ConsolidateRobustness(ingroupAUC_ScGeneRai_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                              minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_GBM_Int <- GetResultRanges(ingroupAUC_ScGeneRai_GBM_Int)
# GBM macrophages
ingroupAUC_ScGeneRai_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirScGeneRai, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirScGeneRai, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_ScGeneRai_GBM_Mac <- ConsolidateRobustness(ingroupAUC_ScGeneRai_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                              minTextDistAsPercentage = 20)
dev.off()
ranges_ScGeneRai_GBM_Mac <- GetResultRanges(ingroupAUC_ScGeneRai_GBM_Mac)
# Consolidate results.
ranges_ScGeneRai <- do.call(rbind,  list(ranges_ScGeneRai_RCC_Mac, ranges_ScGeneRai_RCC_Pod, ranges_ScGeneRai_GBM_Oligo,
                                         ranges_ScGeneRai_GBM_Neur, ranges_ScGeneRai_GBM_Int, ranges_ScGeneRai_GBM_Mac))
ranges_ScGeneRai$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                               rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                               rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_ScGeneRai$Method <- rep("ScGeneRai", 36)
ranges_ScGeneRai <- ranges_ScGeneRai[,c("Method", "CellType", colnames(ranges_ScGeneRai)[1:(ncol(ranges_ScGeneRai)-2)])]


# Read all AUC/ROC results from SINGE.
# RCC Macrophages
ingroupAUC_SINGE_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                                  "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_SINGE_RCC_Mac <- ConsolidateRobustness(ingroupAUC_SINGE_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_RCC_Mac <- GetResultRanges(ingroupAUC_SINGE_RCC_Mac)
# RCC Podocytes
ingroupAUC_SINGE_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_SINGE_RCC_Pod <- ConsolidateRobustness(ingroupAUC_SINGE_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_RCC_Pod <- GetResultRanges(ingroupAUC_SINGE_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_SINGE_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_SINGE_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_SINGE_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                            minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_GBM_Oligo <- GetResultRanges(ingroupAUC_SINGE_GBM_Oligo)
# GBM neurons
ingroupAUC_SINGE_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_SINGE_GBM_Neur <- ConsolidateRobustness(ingroupAUC_SINGE_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_GBM_Neur <- GetResultRanges(ingroupAUC_SINGE_GBM_Neur)
# GBM interneurons
ingroupAUC_SINGE_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_SINGE_GBM_Int <- ConsolidateRobustness(ingroupAUC_SINGE_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_GBM_Int <- GetResultRanges(ingroupAUC_SINGE_GBM_Int)
# GBM macrophages
ingroupAUC_SINGE_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirSINGE, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirSINGE, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_SINGE_GBM_Mac <- ConsolidateRobustness(ingroupAUC_SINGE_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                          minTextDistAsPercentage = 20)
dev.off()
ranges_SINGE_GBM_Mac <- GetResultRanges(ingroupAUC_SINGE_GBM_Mac)
# Consolidate results.
ranges_SINGE <- do.call(rbind,  list(ranges_SINGE_RCC_Mac, ranges_SINGE_RCC_Pod, ranges_SINGE_GBM_Oligo,
                                     ranges_SINGE_GBM_Neur, ranges_SINGE_GBM_Int, ranges_SINGE_GBM_Mac))
ranges_SINGE$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                           rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                           rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_SINGE$Method <- rep("SINGE", 36)
ranges_SINGE <- ranges_SINGE[,c("Method", "CellType", colnames(ranges_SINGE)[1:(ncol(ranges_SINGE)-2)])]


# Read all AUC/ROC results from GRISLI.
rocAucDirGRISLI <- "/Users/tae771/Documents/postdoc/single_cell_CPTAC/GRISLI_results_roc_auc_updatedFERRET/"
# RCC Macrophages
ingroupAUC_GRISLI_RCC_Mac <- lapply(1:length(ingroupNameRCCMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListRCC[i], "_", ingroupNameRCCMac[i], 
                                  "_", outgroupNameRCCMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "RCC_consolidated_macrophage.pdf"))
consolidatedAUCROC_GRISLI_RCC_Mac <- ConsolidateRobustness(ingroupAUC_GRISLI_RCC_Mac, xlab = "Macrophages to Others (RCC)", ylab = "Macrophages to Macrophages (RCC)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_RCC_Mac <- GetResultRanges(ingroupAUC_GRISLI_RCC_Mac)
# RCC Podocytes
ingroupAUC_GRISLI_RCC_Pod <- lapply(1:length(ingroupNameRCCPod), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListRCC[i], "_", ingroupNameRCCPod[i], 
                                  "_", outgroupNameRCCPod[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "RCC_consolidated_Podocyte.pdf"))
consolidatedAUCROC_GRISLI_RCC_Pod <- ConsolidateRobustness(ingroupAUC_GRISLI_RCC_Pod, xlab = "Podocytes to Others (RCC)", ylab = "Podocytes to Podocytes (RCC)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_RCC_Pod <- GetResultRanges(ingroupAUC_GRISLI_RCC_Pod)
# GBM oligodendrocytes
ingroupAUC_GRISLI_GBM_Oligo <- lapply(1:length(ingroupNameGBMOligo), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListGBMOligo[i], "_", ingroupNameGBMOligo[i], 
                                  "_", outgroupNameGBMOligo[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "GBM_consolidated_Oligodendrocyte.pdf"))
consolidatedAUCROC_GRISLI_GBM_Oligo <- ConsolidateRobustness(ingroupAUC_GRISLI_GBM_Oligo, xlab = "Oligodendrocytes to Others (GBM)", ylab = "Oligodendrocytes to Oligodendrocytes (GBM)",
                                                             minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_GBM_Oligo <- GetResultRanges(ingroupAUC_GRISLI_GBM_Oligo)
# GBM neurons
ingroupAUC_GRISLI_GBM_Neur <- lapply(1:length(ingroupNameGBMNeur), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListGBMNeur[i], "_", ingroupNameGBMNeur[i], 
                                  "_", outgroupNameGBMNeur[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "GBM_consolidated_Neuron.pdf"))
consolidatedAUCROC_GRISLI_GBM_Neur <- ConsolidateRobustness(ingroupAUC_GRISLI_GBM_Neur, xlab = "Neurons to Others (GBM)", ylab = "Neurons to Neurons (GBM)",
                                                            minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_GBM_Neur <- GetResultRanges(ingroupAUC_GRISLI_GBM_Neur)
# GBM interneurons
ingroupAUC_GRISLI_GBM_Int <- lapply(1:length(ingroupNameGBMInt), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListGBMInt[i], "_", ingroupNameGBMInt[i], 
                                  "_", outgroupNameGBMInt[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "GBM_consolidated_Interneuron.pdf"))
consolidatedAUCROC_GRISLI_GBM_Int <- ConsolidateRobustness(ingroupAUC_GRISLI_GBM_Int, xlab = "Interneurons to Others (GBM)", ylab = "Interneurons to Interneurons (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_GBM_Int <- GetResultRanges(ingroupAUC_GRISLI_GBM_Int)
# GBM macrophages
ingroupAUC_GRISLI_GBM_Mac <- lapply(1:length(ingroupNameGBMMac), function(i){
  return(ReadRobustnessAUC(paste0(rocAucDirGRISLI, sourceFileListGBMMac[i], "_", ingroupNameGBMMac[i], 
                                  "_", outgroupNameGBMMac[i], "_results.csv")))
})
pdf(paste0(rocAucDirGRISLI, "GBM_consolidated_Macrophage.pdf"))
consolidatedAUCROC_GRISLI_GBM_Mac <- ConsolidateRobustness(ingroupAUC_GRISLI_GBM_Mac, xlab = "Macrophages to Others (GBM)", ylab = "Macrophages to Macrophages (GBM)",
                                                           minTextDistAsPercentage = 20)
dev.off()
ranges_GRISLI_GBM_Mac <- GetResultRanges(ingroupAUC_GRISLI_GBM_Mac)
# Consolidate results.
ranges_GRISLI <- do.call(rbind,  list(ranges_GRISLI_RCC_Mac, ranges_GRISLI_RCC_Pod, ranges_GRISLI_GBM_Oligo,
                                      ranges_GRISLI_GBM_Neur, ranges_GRISLI_GBM_Int, ranges_GRISLI_GBM_Mac))
ranges_GRISLI$CellType <- c(rep("Macrophage (RCC)", 6), rep("Podocyte (RCC)", 6),
                            rep("Oligodendrocyte (GBM)", 6), rep("Neuron (GBM)", 6),
                            rep("Interneuron (GBM)", 6), rep("Macrophage (GBM)", 6))
ranges_GRISLI$Method <- rep("GRISLI", 36)
ranges_GRISLI <- ranges_GRISLI[,c("Method", "CellType", colnames(ranges_GRISLI)[1:(ncol(ranges_GRISLI)-2)])]

# Consolidate all files.
ranges_all <- do.call(rbind, list(ranges_scSGL, ranges_LEAP, ranges_PIDC, ranges_ScGeneRai,
                                  ranges_SINGE, ranges_GRISLI))
write.csv(ranges_all, ranges_all_file)

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