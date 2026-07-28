setwd("/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/netZooR_tmpCommit/netZooR/")
roxygen2::roxygenize()
library("fgsea")

# Initialize directories for each method.
pathway_dir_GRISLI <- NULL
pathway_dir_LEAP <- NULL
pathway_dir_PIDC <- NULL
pathway_dir_scGeneRai <- NULL
pathway_dir_scSGL <- NULL
pathway_dir_SINGE <- NULL
pathway_dir_GRISLI_all <- NULL
pathway_dir_LEAP_all <- NULL
pathway_dir_PIDC_all <- NULL
pathway_dir_scGeneRai_all <- NULL
pathway_dir_scSGL_all <- NULL
pathway_dir_SINGE_all <- NULL
pathway_dir_GRISLI_table <- NULL
pathway_dir_LEAP_table <- NULL
pathway_dir_PIDC_table <- NULL
pathway_dir_scGeneRai_table <- NULL
pathway_dir_scSGL_table <- NULL
pathway_dir_SINGE_table <- NULL

# Load data.
writePathways <- function(pathway_dir, out_dir, table_dir){
  filenameList <- list.files(pathway_dir)
  pathways <- lapply(filenameList, function(filename){
    
    # Get the pathways.
    print(filename)
    pathwayResults <- data.frame(pathway = c(),
                                 ingroup = c(),
                                 outgroup = c(),
                                 sample = c())
    tryCatch({
      pathway <- read.csv(paste0(pathway_dir, filename), row.names = 1)
      sigPathways <- pathway[which(pathway$padj < 0.05), "pathway"]
      
      # Get the cell types.
      underscoreSplit <- strsplit(filename, split = "_")[[1]]
      lastPart <- underscoreSplit[(length(underscoreSplit) - 2):(length(underscoreSplit) - 1)]
      sample <- paste(underscoreSplit[1:(length(underscoreSplit) - 3)], collapse = "_")
      cellTypes <- lastPart[1:2]
      
      # Create data frame.
      pathwayResults <- data.frame(pathway = sigPathways,
                                   ingroup = rep(cellTypes[1], length(sigPathways)),
                                   outgroup = rep(cellTypes[2], length(sigPathways)),
                                   sample = rep(sample, length(sigPathways)))
    }, error = function(cond){
      print(cond)
    })
    
    str(sample)
    str(cellTypes)
    str(pathwayResults)
    return(pathwayResults)
  })
  pathwayDf <- do.call(rbind, pathways)
  write.csv(pathwayDf, out_dir, row.names = FALSE)
  
  cellTypesIngroup <- c(rep("Macrophages", 4), rep("Interneurons", 2), rep("Neurons", 2), rep("Oligodendrocytes", 3), "Podocytes")
  cellTypesOutgroup <- c("Podocytes", "Interneurons", "Neurons", "Oligodendrocytes", 
                         "Macrophages", "Oligodendrocytes", 
                         "Macrophages", "Oligodendrocytes", 
                         "Interneurons", "Neurons", "Macrophages",
                         "Macrophages")
  cellTypeCounts <- c(2, 2, 2, 3, 2, 2, 2, 3, 2, 3, 3, 2)
  for(ii in 1:length(cellTypesIngroup)){
    cellTypeIngroup <- cellTypesIngroup[ii]
    cellTypeOutgroup <- cellTypesOutgroup[ii]
    whichCellType <- intersect(which(pathwayDf$ingroup == cellTypeIngroup), which(pathwayDf$outgroup == cellTypeOutgroup))
    pathwayTable <- table(pathwayDf[whichCellType, "pathway"])
    pathwayTableNorm <- pathwayTable / cellTypeCounts[ii]
    write.csv(pathwayTableNorm, paste0(table_dir, 
                                       cellTypeIngroup, "_", cellTypeOutgroup, ".csv"), row.names = FALSE)
    str(paste(cellTypeIngroup, cellTypeOutgroup))
    write.csv(pathwayTableNorm, paste0(pathway_dir, cellTypeIngroup, "_", cellTypeOutgroup, "_pathways.csv"))
  }
}
writePathways(pathway_dir_GRISLI, pathway_dir_GRISLI_all, pathway_dir_GRISLI_table)
writePathways(pathway_dir_LEAP, pathway_dir_LEAP_all, pathway_dir_LEAP_table)
writePathways(pathway_dir_PIDC, pathway_dir_PIDC_all, pathway_dir_PIDC_table)
writePathways(pathway_dir_scGeneRai, pathway_dir_scGeneRai_all, pathway_dir_scGeneRai_table)
writePathways(pathway_dir_scSGL, pathway_dir_scSGL_all, pathway_dir_scSGL_table)
writePathways(pathway_dir_SINGE, pathway_dir_SINGE_all, pathway_dir_SINGE_table)
