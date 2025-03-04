setwd("/home/ubuntu/netZooR/")
roxygen2::roxygenize()
#library("biomaRt")
library("fgsea")

# Read pathways.
pathways = fgsea::gmtPathways("/home/ubuntu/pathways.gmt")

# Initialize directories for each method.
#pathway_dir <- "/home/ubuntu/SCORPION_pathways"
pathway_dir <- "/home/ubuntu/SCENIC_pathways"
if(!file.exists(pathway_dir)){
  dir.create(pathway_dir)
}

# Load data.
#networks <- LoadResults(resultDirectory = "/home/ubuntu/HTAN_SCORPION/", 
#                        interpretationOfNegative = "poor")
networks <- LoadResults(resultDirectory = "/home/ubuntu/HTAN_SCENIC_processed/", 
                        interpretationOfNegative = "inhibitory", firstColumnIsRowname = FALSE, isTabDelimited = TRUE)

# Helper function to build comparison objects.
foldCount <- 5
getFileNames <- function(rootName, foldCount, ending){
  comparisonList <- unlist(lapply(1:foldCount, function(j){
    retVal <- paste0(paste(rootName, 20, j, sep = "_"), ending)
    print(paste0(networks@directory, "/", retVal))
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

for(i in 1:length(sourceFileList)){
#for(i in 1:2){
  # Compute the AUC and plot the curve.
  #comparisons <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], "logcounts.csv", ingroupName[i], "original.tsv.csv", sep = "_"),
  #                                    ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv.csv"),
  #                                    outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv.csv"),
  #                                    results = networks)
  comparisons <- BuildComparisonObject(sourceNetwork = paste(sourceFileList[i], "logcounts.csv", ingroupName[i], "original.tsv", sep = "_"),
                                      ingroupToCompare = getFileNames(paste(sourceFileList[i], ingroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv"),
                                      outgroupToCompare = getFileNames(paste(sourceFileList[i], outgroupName[i], "logcounts.csv", sep = "_"), 5, ".tsv"),
                                      results = networks)
  
  # Obtain the ingroup and outgroup only and run pathway analysis.
  differentialNetwork <- netZooR::GetDifferentialNetworks(results = networks, comparisons = comparisons)
  colnames(differentialNetwork)[3] <- "Weight"
  str(differentialNetwork)
  differentialNetworkMat <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(differentialNetwork), attr = "Weight"))
  str(differentialNetworkMat)
  differentialPathways <- netZooR::RunPathwayAnalysis(network = differentialNetworkMat,
                                                      pathways = pathways,
  						      geneSymbols = colnames(differentialNetworkMat))
  str(differentialPathways)
  write.csv(differentialPathways, paste0(pathway_dir, "/", paste(sourceFileList[i], 
                                                                 ingroupName[i], outgroupName[i], sep = "_"), "_pathways.csv"))
}
