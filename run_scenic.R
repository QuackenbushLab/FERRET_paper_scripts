# Reference: https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html

# Load SCENIC.
if(!require("SCENIC")){
  devtools::install_github("aertslab/SCENIC") 
}
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
  library(KernSmooth)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(ComplexHeatmap)
})

# Save the SCENIC RMD file (only need to do once).
#vignetteFile <- file.path(system.file('doc', package='SCENIC'), "SCENIC_Running.Rmd")
#file.copy(vignetteFile, "SCENIC_myRun.Rmd")

# Download the files from the legacy site: https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/

# Initialize SCENIC.
cisTarget_database_dir <- NULL
scenicOptions <- initializeScenic(org="hgnc", dbDir=cisTarget_database_dir, nCores=10)


# Load the expression file and set the output directory to the local directory.
args <- commandArgs(trailingOnly = TRUE)
exprMat <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = TRUE))

# Filter the expression network.
dir.create(args[2])
setwd(args[2])
dir.create(paste0(getwd(), "/int"))
dir.create(paste0(getwd(), "/output"))
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

# Run correlation.
runCorrelation(exprMat_filtered, scenicOptions)

# Run GENIE.
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipTsne = TRUE, skipHeatmap = TRUE)

# Binarize activity and save the tSNE plot.
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Save as a LOOM file.
export2loom(scenicOptions, exprMat) # Temporary, to add to loom

# Process the LOOM file to obtain the proper network format.
#loom <- open_loom(getOutputFile("loomFile"))
#regulonAUC <- get_regulons_AUC(loom)
#regulons <- regulonsToGeneLists(get_regulons(loom, column.attr.name="MotifRegulons"))
