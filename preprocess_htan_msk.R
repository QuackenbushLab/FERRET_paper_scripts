# From HTAN, we downloaded Organ Lung NOS, Atlas HTAN MSK, File Format CSV.
# Paper reference: https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00497-9
# We removed the following samples:
# RU 681 (brain metastasis)
# RU 666 (bone metastasis)
# RU 255 (brain metastasis)
# MRM35429319 (spinal fluid ... ?)
# 35379369 (spinal fluid ... ?)
# 1559_1262C P96 (not listed in file)
# RU1141 (not listed in file)
# RU1138 (liver metastasis)

# Required to convert the data to a Bioconductor format
if(!require("SingleCellExperiment")){
  BiocManager::install("SingleCellExperiment")
}
library("SingleCellExperiment")

# Required to get gene symbols and locations.
if(!require("scater")){
  BiocManager::install("scater")
}
library("scater")

# Required to compute size factors.
if(!require("scran")){
  BiocManager::install("scran")
}
library("scran")
# Required to annotate cells.
if(!require("scMRMA")){
  devtools::install_github("JiaLiVUMC/scMRMA")
}
library("scMRMA")

# pip install synapseclient
dataDir <- NULL
# This function is used to find data marked as "NA" when converting to a matrix.

for(file in list.files(dataDir)){
  str(strsplit(file, "_dense.csv")[[1]])
  if(length(strsplit(file, "_dense.csv")[[1]]) == 1){
    data <- t(read.csv(paste0(dataDir, "/", file), header = TRUE, row.names = 1))
    dataMat <- as.matrix(data)
    dataMat <- dataMat[-1,] # Remove cluster annotations.
    
    # Get symbols and locations of genes.
    sce <- SingleCellExperiment(assays = list(counts = dataMat))
    sceLoc <- scater::getBMFeatureAnnos(sce, 
                                        filters = "hgnc_symbol",
                                        attributes = c("hgnc_symbol",
                                                       "start_position", "end_position", "chromosome_name"),
                                        dataset = "hsapiens_gene_ensembl",
                                        host = "useast.ensembl.org")
    print("Finished symbols and locations")
    
    # Get size factors.
    sceSize <- scran::computeSumFactors(sceLoc)
    print("Finshed size factors")
    
    # Normalize.
    sceNorm <- logNormCounts(sceSize)
    print("Finished normalizing")
    
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
    write.csv(qcResults@listData, paste0(dataDir, "/", file, "_QC.csv"))
    
    # Create QC plot.
    png(paste0(dataDir, "/", file, "_features_mito_qc.png"))
    plot(x = qcResults$subsets_mito_sum, y = qcResults$detected,
         main = file, xlab = "Percent Mitochondrial Counts", ylab = "Features By Count")
    dev.off()
    print("Finished QC")
    
    # Save data.
    saveRDS(sceUMAP, paste0(dataDir, "/", file, ".RDS"))
    
    # Filter duplicates.
    sceMatrix <- sceUMAP@assays@data@listData$logcounts
    rownames(sceMatrix) <- rownames(dataMat)
    
    # Assign cell types using scMRMA.
    cellTypes <- scMRMA(input = sceMatrix, species = "Hs", db = "panglaodb", p = 0.05, normalizedData = FALSE, selfDB = NULL, selfClusters = NULL, k = 20)
    sceUMAP@colData@listData$cellTypes <- cellTypes
    saveRDS(sceUMAP, paste0(dataDir, "/", file, "_SCEWithCellTypes.RDS"))
    print(table(cellTypes$multiR$annotationResult$Level3))
    
    # Plot with cell types.
    cellTypes <- as.character(sceUMAP@colData@listData$cellTypes$uniformR$annotationResult$UniformR)
    #allColors <- grDevices::colors()[grep('white', grDevices::colors(), invert = T)]
    #cellTypeColors <- sample(allColors, length(cellTypes))
    #names(cellTypeColors) <- cellTypes
    uniqueCellTypes <- unique(cellTypes)
    colors <- colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(uniqueCellTypes))
    names(colors) <- uniqueCellTypes
    
    # Plot UMAP.
    umap <- sceUMAP@int_colData@listData$reducedDims@listData$UMAP
    sceDF <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2], cellType = cellTypes, color = colors[cellTypes])
    offset <- max(sceDF$UMAP1) - min(sceDF$UMAP1) / 3
    title <-  strsplit(file, ".", fixed = TRUE)[[1]][1]
    png(paste0(dataDir, "/", file, "_cellTypesUMAP.png"))
    par(mar=c(4,4,4,15), xpd = TRUE)
    plot(sceDF$UMAP1, sceDF$UMAP2, pch = '.', col = sceDF$color, xlab = "UMAP1",
         ylab = "UMAP2", main = title)
    legend("topright", inset = c(-0.8, 0), legend=uniqueCellTypes, fill=colors, title="Cell Type")
    dev.off()
    
    # Plot cell types with respect to QC.
    png(paste0(dataDir, "/", file, "_cellTypesQC.png"))
    par(mar=c(4,4,4,15), xpd = TRUE)
    plot(qcResults$subsets_mito_sum, qcResults$detected,
         main = title, col = sceDF$color, pch = '.', xlab = "Percent Mitochondrial Counts", ylab = "Features By Count",
         xlim = c(0, 80))
    legend("topright", inset = c(-0.7, 0), legend=uniqueCellTypes, fill=colors, title="Cell Type")
    dev.off()
  }
}

# Do QC filtering and write out the cell type specific files.
for(file in list.files(dataDir)){
  if(length(strsplit(file, "_dense.csv")[[1]]) == 1 && strsplit(file, "_dense.csv")[[1]][1] != file){
    sce <- readRDS(paste0(dataDir,"/", file, "_SCEWithCellTypes.RDS"))
    qcResults <- read.csv(paste0(dataDir, "/", file, "_QC.csv"))
    
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
    saveRDS(sceFeat, paste0(dataDir, "/", file, "_SCEQC.RDS"))
    
    # Create a directory for all cell-type inputs.
    cellTypeInputDir <- paste0(dataDir, "/cellTypeSpecificMatrices")
    if(!file.exists(cellTypeInputDir)){
      dir.create(cellTypeInputDir)
    }
    
    # For each cell type, write out the file. Write out the original file too.
    cellTypes <- unique(sceFeat@colData@listData$cellTypes$uniformR$annotationResult)
    for(cellType in cellTypes){
      write.csv(sceFeat@assays@data@listData$logcounts[,which(sceFeat@colData@listData$cellTypes$uniformR$annotationResult == cellType)],
                paste0(cellTypeInputDir, "/", file, "_", cellType, "_logcounts.csv"))
    }
    write.csv(sceFeat@assays@data@listData$logcounts,
              paste0(cellTypeInputDir, "/", file, "_logcounts.csv"))
  }
}
