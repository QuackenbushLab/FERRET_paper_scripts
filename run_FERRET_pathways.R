library("FERRET")

# Variables
scGeneRai_results <- NULL
GRISLI_results <- NULL
scSGL_results <- NULL
PIDC_results <- NULL
LEAP_results <- NULL
SINGE_results <- NULL
pathway_dir <- NULL
GMT <- NULL

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

# Read PCs
pca <- lapply(uniqueSamples, function(sample){
  file <- readRDS(paste0("/home/ubuntu/TCGA_single_cell_pca_results/",
                         sample, "_pca.RDS"))
  rownames(file$rotation) <- names(file$center)
  return(file$rotation)
})
names(pca) <- uniqueSamples

# Map all genes to symbol names.
geneNames <- unique(unlist(lapply(pca, function(pc){
  return(rownames(pc))
})))
# Need to install an old version of dbplyr to get this to work.
#install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)
symbolMapping <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol", "ensembl_gene_id"),
                       values= geneNames,mart= mart)
symbolMappingUniqueEnsembl <- unique(symbolMapping$ensembl_gene_id)
symbolMappingUniqueEnsemblGene <- unlist(lapply(symbolMappingUniqueEnsembl, function(ensembl){
  return(symbolMapping[which(symbolMapping$ensembl_gene_id == ensembl), "hgnc_symbol"][1])
}))
symbolMappingUnique <- data.frame(ensembl = symbolMappingUniqueEnsembl,
                                  gene = symbolMappingUniqueEnsemblGene)
rownames(symbolMappingUnique) <- symbolMappingUnique$ensembl

# Read pathways (downloaded from MSIGDB).
pathways = gmtPathways(GMT)

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
      
      differentialNetworkScSGL <- netZooR::GetDifferentialNetworks(results = networksScSGL, comparisons = comparisonsScSGL)
      differentialNetworkGenesScSGL <- netZooR::GetFullNetworkFromPCANetwork(pcNetwork = differentialNetworkScSGL, PC = pca[[sourceFileList[i]]])
      differentialNetworkGenesSymbolScSGL <- differentialNetworkGenesScSGL
      differentialNetworkGenesSymbolScSGL[,2] <- unlist(lapply(differentialNetworkGenesScSGL[,2], function(gene){
        return(symbolMappingUnique[gene,"gene"])
      }))
      differentialPathwaysScSGL <- netZooR::RunPathwayAnalysis(network = differentialNetworkGenesSymbolScSGL,
                                                          pathways = pathways)
      write.csv(differentialPathwaysScSGL, paste0(pathway_dir, "/", paste(sourceFileList[i], 
                                                                     ingroupName[i], outgroupName[i], sep = "_"), "_pathways.csv"))
    }, error = function(cond){
      print(paste("Skipping iteration", i))
    })
}