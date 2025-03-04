networkDir <- NULL
resultDir <- NULL
pathwayFile <- NULL
differentialPathwaysDir <- NULL

# Load the networks. (Note that we copied them.)
networks <- LoadResults(networkDir, firstColumnIsRowname = FALSE)

# Build the in-group (B-lymphocyte) and out-group (fibroblast) objects.
comparisonsBF <- BuildComparisonObject(sourceNetwork = "B_lymphocytes_networkAll.csv",
                                       ingroupToCompare = paste0("B_lymphocytes_networkWithout", 1:6, ".csv"),
                                       outgroupToCompare = paste0("fibroblasts_networkWithout", 1:2, ".csv"),
                                       results = networks)

# Run FERRET.
pdf(paste0(resultDir, "/RAUC/rocBF.pdf"))
result <- ComputeRobustnessAUC(results = networks,comparisons = comparisonsBF,xlab = "B-lymphocytes to Fibroblasts",
                               ylab = "B-lymphocytes to B-lymphocytes")
dev.off()

# Build the in-group (fibroblast) and out-group (B-lymphocyte) objects.
comparisonsBF <- BuildComparisonObject(sourceNetwork = "fibroblasts_networkAll.csv",
                                       ingroupToCompare = paste0("fibroblasts_networkWithout", 1:2, ".csv"),
                                       outgroupToCompare = paste0("B_lymphocytes_networkWithout", 1:6, ".csv"),
                                       results = networks)

# Run FERRET.
pdf(paste0(resultDir, "/RAUC/rocFB.pdf"))
result <- ComputeRobustnessAUC(results = networks,comparisons = comparisonsBF,xlab = "Fibroblasts to B-lymphocytes",
                               ylab = "Fibroblasts to Fibroblasts")
dev.off()

# Do pathway analysis.
library("fgsea")
pathways = gmtPathways(pathwayFile)
comparisonsBF <- BuildComparisonObject(sourceNetwork = "B_lymphocytes_networkAll.csv",
                                       ingroupToCompare = paste0("B_lymphocytes_networkWithout", 1:6, ".csv"),
                                       outgroupToCompare = paste0("fibroblasts_networkWithout", 1:2, ".csv"),
                                       results = networks)
differentialNetwork <- netZooR::GetDifferentialNetworks(results = networks, comparisons = comparisonsBF)
colnames(differentialNetwork)[3] <- "weight"
differentialNetworkAdj <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(differentialNetwork), attr = "weight")
differentialPathways <- netZooR::RunPathwayAnalysis(network = differentialNetworkAdj,
                                                    pathways = pathways, geneSymbols = colnames(differentialNetworkAdj))
write.csv(differentialPathways, paste0(differentialPathwaysDir, "/B_lymphocyte_pathways/pw.csv"))
comparisonsFB <- BuildComparisonObject(sourceNetwork = "fibroblasts_networkAll.csv",
                                       ingroupToCompare = paste0("fibroblasts_networkWithout", 1:2, ".csv"),
                                       outgroupToCompare = paste0("B_lymphocytes_networkWithout", 1:6, ".csv"),
                                       results = networks)
differentialNetwork <- netZooR::GetDifferentialNetworks(results = networks, comparisons = comparisonsFB)
colnames(differentialNetwork)[3] <- "weight"
differentialNetworkAdj <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(differentialNetwork), attr = "weight")
differentialPathways <- netZooR::RunPathwayAnalysis(network = differentialNetworkAdj,
                                                    pathways = pathways, geneSymbols = colnames(differentialNetworkAdj))
write.csv(differentialPathways, paste0(differentialPathwaysDir, "/fibroblast_pathways/pw.csv"))
