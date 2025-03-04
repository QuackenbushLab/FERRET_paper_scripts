# Reference: https://bioconductor.org/books/3.13/OSCA.basic/
# Reference 2: https://kieranrcampbell.github.io/r-workshop-march-2019/
# In this version, to perfrom the cell annotation, we use scMRMA: https://github.com/JiaLiVUMC/scMRMA
if(!require("devtools")){
  install.packages("devtools")
}
# Required to read in the data
if(!require("loomR")){
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop") 
}

# Required to convert the data to a Bioconductor format
if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
}

# Required to assign cell type
if(!require("SingleR")){
  BiocManager::install("SingleR")
}
if(!require("scater")){
  BiocManager::install("scater")
}
if(!require("scran")){
  BiocManager::install("scran")
}
if(!require("scMRMA")){
  devtools::install_github("JiaLiVUMC/scMRMA")
}
if(!require("RColorBrewer")){
  install.packages("RColorBrewer")
}
if(!require("dichromat")){
  install.packages("dichromat")
}

# Variables to set.
pythonpath <- NULL
dataDir <- NULL
biospecimenDir <- NULL
clinicalDir <- NULL

# Required to plot cell type.
#if(!require("SeuratObject")){
# install.packages("SeuratObject")
#}
#if(!require("Seurat")){
#  remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
#}
Sys.setenv(RETICULATE_PYTHON=)
library("tidyverse")
library("pheatmap")
library("DT")
library("devtools")
library("hdf5r")
library("loomR")
library("SingleR")
library("scran")
library("SingleCellExperiment")
library(tensorflow)
library("readr")
library("scMRMA")
library("RColorBrewer")
library("dichromat")
library("data.table")
#library("SeuratObject")
#library("Seurat")

# Read in all sample and clinical data and reference databases.
sampleSheet <- fread(paste0(dataDir, "/", biospecimenDir, "/sample.tsv"), sep = "\t",
                          header = TRUE)
uniqueSamples <- unique(sampleSheet$sample_submitter_id)
clinical <- fread(paste0(dataDir, "/", clinicalDir, "/clinical.tsv"), sep = "\t",
                       header = TRUE)

# For each sample, perform preprocessing and QC.
for(i in 1:length(uniqueSamples)){
  sample <- uniqueSamples[i]
  print(paste("Starting processing for sample", sample))
  # Separate out the files.
  filesForSample <- paste(dataDir, paste(sampleSheet[which(sampleSheet$Sample.ID == sample), "File.ID"],
                          sampleSheet[which(sampleSheet$Sample.ID == sample), "File.Name"], sep = "."), sep = "/gene_expression/")
  
  # Convert Loom to SingleCellExperiment for compatibility with the package.
  loomFile <- filesForSample[which(grepl("loom", filesForSample, fixed = TRUE) == TRUE)]
  loom <- loomR::connect(filename = loomFile, skip.validate = TRUE)
  loomMatrix <- loom[["matrix"]][,]
  colnames(loomMatrix) <- unlist(lapply(loom[["row_attrs/Gene"]][], function(gene){
    return(strsplit(gene, ".", fixed = TRUE)[[1]][1])
  }))
  rownames(loomMatrix) <- loom[["col_attrs/CellID"]][]
  sce <- SingleCellExperiment(assays = list(counts = t(loomMatrix)))

  # Get symbols and locations of genes.
  sceLoc <- scater::getBMFeatureAnnos(sce, 
                    filters = "ensembl_gene_id",
                    attributes = c("ensembl_gene_id", "hgnc_symbol",
                                   "start_position", "end_position", "chromosome_name"),
                    dataset = "hsapiens_gene_ensembl")
  print("Finished symbols and locations")
  
  # Get size factors.
  sceSize <- scran::computeSumFactors(sceLoc)
  print("Finshed size factors")
  
  # Normalize.
  sceNorm <- logNormCounts(sceSize)
  print("Finished normalizing")
  
  # Save the histogram.
  png(paste0(dataDir, "/", sample, "_histogram.png"))
  hist(sceNorm@assays@data@listData$counts, main = sample,
       xlab = "Gene Expression")
  dev.off()
  print("Finished histogram")
  
  # Run dimensionality reduction.
  scePCA <- scater::runPCA(sceNorm)
  sceUMAP <- scater::runUMAP(scePCA)
  print("Finished dimensionality reduction")
  
  # Calculate and save the QC metrics, calculated using
  # mitochondrial and ribosomal genes.
  mt_genes <- which(rowData(sceUMAP)$chromosome_name == "MT")
  ribo_genes <- grepl("^RP[LS]", rowData(sceUMAP)$Symbol)
  feature_ctrls <- list(mito = rownames(sceUMAP)[mt_genes],
                        ribo = rownames(sceUMAP)[ribo_genes])
  qcResults <- scater::perCellQCMetrics(sceUMAP, subsets = feature_ctrls)
  write.csv(qcResults@listData, paste0(dataDir, "/", sample, "_QC.csv"))
  
  # Create QC plot.
  png(paste0(dataDir, "/", sample, "_features_mito_qc.png"))
  plot(x = qcResults$subsets_mito_sum, y = qcResults$detected,
       main = sample, xlab = "Percent Mitochondrial Counts", ylab = "Features By Count")
  dev.off()
  print("Finished QC")
  
  # Save data.
  saveRDS(sceUMAP, paste0(dataDir, "/", sample, ".RDS"))
}


# For each sample, obtain cell type.
for(i in 1:length(uniqueSamples)){
  sample <- uniqueSamples[i]
  print(paste("Obtaining cell type for sample", sample))
  
  # Read data and format as matrix.
  sce <- readRDS(paste0(dataDir, "/", sample, ".RDS"))
  sceMatrix <- sce@assays@data@listData$counts
  rownames(sceMatrix) <- sce@rowRanges@elementMetadata@listData$hgnc_symbol
  
  # Filter duplicates.
  whichHasSymbol <- which(sce@rowRanges@elementMetadata@listData$hgnc_symbol != "")
  sceMatrixNoBlanks <- sceMatrix[whichHasSymbol,]
  sceWhichHasSymbol <- sce
  sceWhichHasSymbol@int_elementMetadata@nrows <- length(whichHasSymbol)
  sceWhichHasSymbol@int_elementMetadata@listData$rowPairs@nrows <- length(whichHasSymbol)
  sceWhichHasSymbol@rowRanges@elementMetadata@nrows <- length(whichHasSymbol)
  sceWhichHasSymbol@rowRanges@elementMetadata@listData$ensembl_gene_id <- sce@rowRanges@elementMetadata@listData$ensembl_gene_id[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@elementMetadata@listData$hgnc_symbol <- sce@rowRanges@elementMetadata@listData$hgnc_symbol[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@elementMetadata@listData$start_position <- sce@rowRanges@elementMetadata@listData$start_position[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@elementMetadata@listData$end_position <- sce@rowRanges@elementMetadata@listData$end_position[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@elementMetadata@listData$chromosome_name <- sce@rowRanges@elementMetadata@listData$chromosome_name[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@partitioning@end <- sce@rowRanges@partitioning@end[whichHasSymbol]
  sceWhichHasSymbol@rowRanges@partitioning@NAMES <- sce@rowRanges@partitioning@NAMES[whichHasSymbol]
  sceWhichHasSymbol@assays@data@listData$counts <- sce@assays@data@listData$counts[whichHasSymbol,]
  sceWhichHasSymbol@assays@data@listData$logcounts <- sce@assays@data@listData$logcounts[whichHasSymbol,]
  sceWhichHasSymbol@elementMetadata@nrows <- length(whichHasSymbol)

  # Filter duplicates.
  rownameTable <- table(rownames(sceMatrixNoBlanks))
  dupIdxToRemove <- lapply(names(rownameTable)[which(rownameTable > 1)], function(gene){
    whichGene <- which(rownames(sceMatrixNoBlanks) == gene)
    sumByGene <- rowSums(sceMatrixNoBlanks[whichGene,])
    lowCoverageIdx <- which.min(sumByGene)
    return(whichGene[lowCoverageIdx])
  })
  dedupToKeep <- setdiff(1:nrow(sceMatrixNoBlanks), unlist(dupIdxToRemove))
  sceMatrixDedup <- sceMatrixNoBlanks[dedupToKeep,]
  sceWhichDedup <- sceWhichHasSymbol
  sceWhichDedup@int_elementMetadata@nrows <- length(dedupToKeep)
  sceWhichDedup@int_elementMetadata@listData$rowPairs@nrows <- length(dedupToKeep)
  sceWhichDedup@rowRanges@elementMetadata@nrows <- length(dedupToKeep)
  sceWhichDedup@rowRanges@elementMetadata@listData$ensembl_gene_id <- sceWhichHasSymbol@rowRanges@elementMetadata@listData$ensembl_gene_id[dedupToKeep]
  sceWhichDedup@rowRanges@elementMetadata@listData$hgnc_symbol <- sceWhichHasSymbol@rowRanges@elementMetadata@listData$hgnc_symbol[dedupToKeep]
  sceWhichDedup@rowRanges@elementMetadata@listData$start_position <- sceWhichHasSymbol@rowRanges@elementMetadata@listData$start_position[dedupToKeep]
  sceWhichDedup@rowRanges@elementMetadata@listData$end_position <- sceWhichHasSymbol@rowRanges@elementMetadata@listData$end_position[dedupToKeep]
  sceWhichDedup@rowRanges@elementMetadata@listData$chromosome_name <- sceWhichHasSymbol@rowRanges@elementMetadata@listData$chromosome_name[dedupToKeep]
  sceWhichDedup@rowRanges@partitioning@end <- sceWhichHasSymbol@rowRanges@partitioning@end[dedupToKeep]
  sceWhichDedup@rowRanges@partitioning@NAMES <- sceWhichHasSymbol@rowRanges@partitioning@NAMES[dedupToKeep]
  sceWhichDedup@assays@data@listData$counts <- sceWhichHasSymbol@assays@data@listData$counts[dedupToKeep,]
  sceWhichDedup@assays@data@listData$logcounts <- sceWhichHasSymbol@assays@data@listData$logcounts[dedupToKeep,]
  sceWhichDedup@elementMetadata@nrows <- length(dedupToKeep)

  # Assign cell types using scMRMA.
  cellTypes <- scMRMA(input = sceMatrixDedup, species = "Hs", db = "panglaodb", p = 0.05, normalizedData = FALSE, selfDB = NULL, selfClusters = NULL, k = 20)
  sceWhichDedup@colData@listData$cellTypes <- cellTypes
  saveRDS(sceWhichDedup, paste0(dataDir, "/", sample, "_SCEWithCellTypes.RDS"))
  print(table(cellTypes$multiR$annotationResult$Level3))
}

# Map each cell type to a color and plot.
brainList <- c("Temporal lobe", "Brain, NOS", "Occipital lobe", "Frontal lobe")
cellTypesBrain <- unique(unlist(lapply(1:length(uniqueSamples), function(i){
  sample <- uniqueSamples[i]
  retval <- c()
  if(clinical[which(clinical$case_submitter_id == paste(strsplit(sample, "-", fixed = TRUE)[[1]][1:2], collapse = "-")), 
                    "tissue_or_organ_of_origin"] %in% brainList){
    print(sample)
    sce <- readRDS(paste0(dataDir, "/", sample, "_SCEWithCellTypes.RDS"))
    cellTypes <- as.character(sce@colData@listData$cellTypes$uniformR$annotationResult$UniformR)
    retval <- cellTypes
    print(unique(retval))
  }
  return(retval)
})))
cellTypesKidney <- unique(unlist(lapply(1:length(uniqueSamples), function(i){
  sample <- uniqueSamples[i]
  retval <- c()
  if(clinical[which(clinical$case_submitter_id == paste(strsplit(sample, "-", fixed = TRUE)[[1]][1:2], collapse = "-")), 
              "tissue_or_organ_of_origin"] == "Kidney, NOS"){
    print(sample)
    sce <- readRDS(paste0(dataDir, "/", sample, "_SCEWithCellTypes.RDS"))
    cellTypes <- as.character(sce@colData@listData$cellTypes$uniformR$annotationResult$UniformR)
    retval <- cellTypes
    print(unique(retval))
  }
  return(retval)
})))
cellTypesAll <- unique(c(cellTypesBrain, cellTypesKidney))
# Solution to generate many colors from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
allColors <- grDevices::colors()[grep('white', grDevices::colors(), invert = T)]
cellTypeColors <- sample(allColors, length(cellTypesAll))
names(cellTypeColors) <- cellTypesAll
for(i in 1:length(uniqueSamples)){
  # Retrieve data.
  sample <- uniqueSamples[i]
  sce <- readRDS(paste0(dataDir, "/preprocessed_data_files/", sample, "_SCEWithCellTypes.RDS"))
  umap <- sce@int_colData@listData$reducedDims@listData$UMAP
  cellTypes <- as.character(sce@colData@listData$cellTypes$uniformR$annotationResult$UniformR)
  uniqueCellTypes <- unique(cellTypes)
  colors <- colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(uniqueCellTypes))
  names(colors) <- uniqueCellTypes
  
  # Plot UMAP.
  sceDF <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2], cellType = cellTypes, color = colors[cellTypes])
  #offset <- max(sceDF$UMAP1) - min(sceDF$UMAP1) / 3
  title <-  paste(sample, 
                  clinical[which(clinical$case_submitter_id == paste(strsplit(sample, "-", 
                                                                              fixed = TRUE)[[1]][1:2], 
                                                                     collapse = "-")), 
                           "primary_diagnosis"])
  #png(paste0(dataDir, "/", sample, "_cellTypesUMAP.png"))
  #par(mar=c(4,4,4,15), xpd = TRUE)
  #plot(sceDF$UMAP1, sceDF$UMAP2, pch = '.', col = sceDF$color, xlab = "UMAP1",
  #     ylab = "UMAP2", main = title)
  #legend("topright", inset = c(-0.8, 0), legend=uniqueCellTypes, fill=colors, title="Cell Type")
  #dev.off()
  
  # Plot cell types with respect to QC.
  qcResults <- read.csv(paste0(dataDir, "/preprocessed_data_files/", sample, "_QC.csv"))
  str(qcResults)
  png(paste0(dataDir, "/cell_type_plots_qc/", sample, "_cellTypesQC.png"))
  par(mar=c(4,4,4,15), xpd = TRUE)
  plot(qcResults$subsets_mito_sum, qcResults$detected,
       main = title, col = sceDF$color, pch = '.', xlab = "Percent Mitochondrial Counts", ylab = "Features By Count",
       xlim = c(0, 80))
  legend("topright", inset = c(-0.7, 0), legend=uniqueCellTypes, fill=colors, title="Cell Type")
  dev.off()
}

# Output cell-type-specific, sample-specific, data sets.
uniqueSamples =  c("C3N-01814-01", "C3N-01904-02", "C3N-02181-02")
for(i in 1:length(uniqueSamples)){
  sample <- uniqueSamples[i]
  #tryCatch({
    sce <- readRDS(paste0(dataDir, "/preprocessed_data_files/", sample, "_SCEWithCellTypes.RDS"))
    qcResults <- read.csv(paste0(dataDir, "/preprocessed_data_files/", sample, "_QC.csv"))
    
    # Remove cells below a features threshold
    # and above a mitochondrial percent threshold.
    sceFeat <- sce
    whichFeatureAndMitoCutoff <- intersect(which(qcResults$detected > 1000),
                                           which(qcResults$subsets_mito_percent < 50))
    sceFeat@assays@data@listData$counts <- sce@assays@data@listData$counts[,whichFeatureAndMitoCutoff]
    sceFeat@assays@data@listData$logcounts <- sce@assays@data@listData$logcounts[,whichFeatureAndMitoCutoff]
    sceFeat@colData@nrows <- length(whichFeatureAndMitoCutoff)
    sceFeat@colData@rownames <- sce@colData@rownames[whichFeatureAndMitoCutoff]
    sceFeat@colData@listData$sizeFactor <- sce@colData@listData$sizeFactor[whichFeatureAndMitoCutoff]
    sceFeat@colData@listData$cellTypes$multiR$annotationResult <- sce@colData@listData$cellTypes$multiR$annotationResult[whichFeatureAndMitoCutoff,]
    sceFeat@colData@listData$cellTypes$multiR$meta <- sce@colData@listData$cellTypes$multiR$meta[whichFeatureAndMitoCutoff,]
    sceFeat@colData@listData$cellTypes$uniformR$annotationResult <- sce@colData@listData$cellTypes$uniformR$annotationResult[whichFeatureAndMitoCutoff,]
    sceFeat@colData@listData$cellTypes$uniformR$meta <- sce@colData@listData$cellTypes$uniformR$meta[whichFeatureAndMitoCutoff,]
    sceFeat@int_colData@nrows <- length(whichFeatureAndMitoCutoff)
    sceFeat@int_colData@rownames <- sce@int_colData@rownames[whichFeatureAndMitoCutoff]
    sceFeat@int_colData@listData$sizeFactor <- sce@int_colData@listData$sizeFactor[whichFeatureAndMitoCutoff]
    sceFeat@int_colData@listData$reducedDims@nrows <- length(whichFeatureAndMitoCutoff)
    sceFeat@int_colData@listData$reducedDims@listData$PCA <- sce@int_colData@listData$reducedDims@listData$PCA[whichFeatureAndMitoCutoff,]
    sceFeat@int_colData@listData$reducedDims@listData$UMAP <- sce@int_colData@listData$reducedDims@listData$UMAP[whichFeatureAndMitoCutoff,]
    sceFeat@int_colData@listData$altExps@nrows <- length(whichFeatureAndMitoCutoff)
    sceFeat@int_colData@listData$colPairs@nrows <- length(whichFeatureAndMitoCutoff)
    sceFeat@colData@listData$qc <- qcResults[whichFeatureAndMitoCutoff,]
    #saveRDS(sceFeat, paste0(dataDir, "/preprocessed_data_files/", sample, "_SCEQC.RDS"))
    
    # For brain cells, keep neurons, interneurons, oligodendrocytes, microglia,
    # and T cells.
    #if(clinical[which(clinical$case_submitter_id == paste(strsplit(sample, "-", fixed = TRUE)[[1]][1:2], collapse = "-")), 
    #            "tissue_or_organ_of_origin"] %in% brainList){
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Neurons")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Neurons_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Neurons")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Neurons_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Interneurons")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Interneurons_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Interneurons")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Interneurons_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Oligodendrocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Oligodendrocytes_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Oligodendrocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Oligodendrocytes_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Macrophages")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Macrophages_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Macrophages")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Macrophages_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Microglia")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Microglia_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Microglia")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Microglia_logcounts.csv"))
      write.csv(sceFeat@assays@data@listData$counts,
                paste0(dataDir, "/preprocessed_data_files/", sample, "_counts.csv"))
    #}
    
    # For kidney cells, keep podocytes, proximal tubule cells, T cells, and macrophages.
    #if(clinical[which(clinical$case_submitter_id == paste(strsplit(sample, "-", fixed = TRUE)[[1]][1:2], collapse = "-")), 
    #            "tissue_or_organ_of_origin"] == "Kidney, NOS"){
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Podocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Podocytes_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Podocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Podocytes_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Proximaltubulecells")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Proximaltubulecells_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Proximaltubulecells")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Proximaltubulecells_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Oligodendrocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Oligodendrocytes_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Oligodendrocytes")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Oligodendrocytes_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Macrophages")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Macrophages_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Macrophages")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Macrophages_logcounts.csv"))
      # write.csv(sceFeat@assays@data@listData$counts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Tcells")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Tcells_counts.csv"))
      # write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == "Tcells")],
      #           paste0(dataDir, "/preprocessed_data_files/", sample, "_Tcells_logcounts.csv"))
      #write.csv(sceFeat@assays@data@listData$counts,
    #            paste0(dataDir, "/preprocessed_data_files/", sample, "_counts.csv"))
    #}
  #}, error = function(cond){print(cond)})
}

