library("FERRET")
library("fgsea")

# Initialize directories for each method.
pathway_dir_SCORPION <- NULL
pathway_dir_SCENIC <- NULL
pathway_dir_SCORPION_all <- NULL
pathway_dir_SCENIC_all <- NULL
pathway_dir_SCORPION_table <- NULL
pathway_dir_SCENIC_table <- NULL

# Load data.
filenameList_SCORPION <- list.files(pathway_dir_SCORPION)
pathways_SCORPION <- lapply(filenameList_SCORPION, function(filename){
  
  # Get the pathways.
  pathway <- read.csv(paste0(pathway_dir, filename), row.names = 1)
  sigPathways <- pathway[which(pathway$padj < 0.05), "pathway"]
  
  # Get the cell types.
  underscoreSplit <- strsplit(filename, split = "_")[[1]]
  lastPart <- underscoreSplit[(length(underscoreSplit) - 4):(length(underscoreSplit) - 1)]
  sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 3)], collapse = "_")
  cellTypes <- lastPart[3:4]
  if(lastPart[2] == "cells"){
    cellTypes[1] <- paste(lastPart[1], lastPart[2])
    sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 5)], collapse = "_")
  }
  if(lastPart[3] == "cells"){
    cellTypes[1] <- paste(lastPart[2], lastPart[3])
    sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 4)], collapse = "_")
  }
  if(lastPart[4] == "cells"){
    if(lastPart[2] != "cells"){
      cellTypes[1] <- lastPart[2] 
      sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 4)], collapse = "_")
    }else{
      sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 5)], collapse = "_")
    }
    cellTypes[2] <- paste(lastPart[3], lastPart[4])
  }
  
  # Create data frame.
  pathwayResults <- data.frame(pathway = sigPathways,
                               ingroup = rep(cellTypes[1], length(sigPathways)),
                               outgroup = rep(cellTypes[2], length(sigPathways)),
                               sample = rep(sample, length(sigPathways)))
  return(pathwayResults)
})
pathwayDf_SCORPION <- do.call(rbind, pathways_SCORPION)
filenameList_SCENIC <- list.files(pathway_dir_SCENIC)
pathways_SCENIC <- lapply(filenameList_SCENIC, function(filename){
  
  # Get the pathways.
  pathway <- read.csv(paste0(pathway_dir, filename), row.names = 1)
  sigPathways <- pathway[which(pathway$padj < 0.05), "pathway"]
  
  # Get the cell types.
  underscoreSplit <- strsplit(filename, split = "_")[[1]]
  lastPart <- underscoreSplit[(length(underscoreSplit) - 4):(length(underscoreSplit) - 1)]
  sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 3)], collapse = "_")
  cellTypes <- lastPart[3:4]
  if(lastPart[2] == "cells"){
    cellTypes[1] <- paste(lastPart[1], lastPart[2])
    sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 5)], collapse = "_")
  }
  if(lastPart[3] == "cells"){
    cellTypes[1] <- paste(lastPart[2], lastPart[3])
    sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 4)], collapse = "_")
  }
  if(lastPart[4] == "cells"){
    if(lastPart[2] != "cells"){
      cellTypes[1] <- lastPart[2] 
      sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 4)], collapse = "_")
    }else{
      sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 5)], collapse = "_")
    }
    cellTypes[2] <- paste(lastPart[3], lastPart[4])
  }
  
  # Create data frame.
  pathwayResults <- data.frame(pathway = sigPathways,
                               ingroup = rep(cellTypes[1], length(sigPathways)),
                               outgroup = rep(cellTypes[2], length(sigPathways)),
                               sample = rep(sample, length(sigPathways)))
  return(pathwayResults)
})
pathwayDf_SCENIC <- do.call(rbind, pathways_SCENIC)

write.csv(pathwayDf_SCORPION, pathway_dir_SCORPION_all, row.names = FALSE)
write.csv(pathwayDf, pathway_dir_SCENIC_all, row.names = FALSE)

#Break it down by cell type.
breakDown <- function(pathwayDf, pathwayTable)
{
  comparisonTypes <- paste(pathwayDf$ingroup, pathwayDf$outgroup, sep = "_")
  cellTypes <- unique(comparisonTypes)
  for(cellType in cellTypes){
    whichCellType <- which(comparisonTypes == cellType)
    pathwayTable <- table(pathwayDf[whichCellType, "pathway"])
    comparisonCount <- length(unique(paste(comparisonTypes[whichCellType], pathwayDf[whichCellType, "sample"])))
    pathwayTableNorm <- pathwayTable / comparisonCount
    write.csv(pathwayTableNorm, paste0(pathwayTable,
                                       cellType, ".csv"), row.names = FALSE)
    str(cellType)
    str(comparisonCount)
    str(pathwayTableNorm)
  }
  cellTypesIngroup <- c(rep("Endothelial cells", 4), rep("Macrophages", 5), rep("B cells", 3), rep("NK cells", 3), rep("Fibroblasts", 4),
                        rep("T cells", 5))
  cellTypesOutgroup <- c("Macrophages", "NK cells", "Fibroblasts", "T cells", 
                         "Endothelial cells", "NK cells", "B cells", "Fibroblasts", "T cells", 
                         "Macrophages", "Fibroblasts", "T cells",
                         "Endothelial cells", "Macrophages", "T cells",
                         "Endothelial cells", "B cells", "Macrophages", "T cells",
                         "Endothelial cells", "B cells", "Macrophages", "NK cells", "Fibroblasts")
  cellTypeCounts <- c(4, 1, 1, 1, 4, 2, 1, 3, 5, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 5, 1, 2)
  for(ii in 1:length(cellTypesIngroup)){
    cellTypeIngroup <- cellTypesIngroup[ii]
    cellTypeOutgroup <- cellTypesOutgroup[ii]
    whichCellType <- intersect(which(pathwayDf$ingroup == cellTypeIngroup), which(pathwayDf$outgroup == cellTypeOutgroup))
    pathwayTable <- table(pathwayDf[whichCellType, "pathway"])
    pathwayTableNorm <- pathwayTable / cellTypeCounts[ii]
    str(paste(cellTypeIngroup, cellTypeOutgroup))
    write.csv(pathwayTableNorm, paste0(pathway_dir, cellTypeIngroup, "_", cellTypeOutgroup, "_pathways.csv"))
  }
}
breakDown(pathwayDf_SCORPION, pathway_dir_SCORPION_table)
breakDown(pathwayDf_SCENIC, pathway_dir_SCENIC_table)
