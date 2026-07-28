# Load SCORPION.
if(!require("SCORPION")){
  install.packages("SCORPION")
}
library("SCORPION")
if(!require("Seurat")){
  install.packages("Seurat")
}
library("Seurat")
if(!require(reshape2)){
  install.packages(reshape2)
}
library(reshape2)

# Variables
motif_file <- NULL
ppi_file <- NULL

# Get parameters.
args <- commandArgs(trailingOnly = TRUE)
print(paste0(args[2], ".csv"))
# Load PANDA motif and PPI priors.
motifs <- read.table(motif_files,sep = "\t")
ppi <- read.table(ppi_file,sep = "\t")
colnames(ppi) <- c("protein1", "protein2", "combined_score")

# Load file.
expression <- as.matrix(read.table(args[1], sep = "\t", row.names = 1, header = TRUE))
str(expression)
# Prepare for SCORPION.
scorpionInput <- list(tf = motifs, ppi = ppi, gex = expression)

# Run SCORPION.
scorpionOutput <- scorpion(tfMotifs = scorpionInput$tf,
                           gexMatrix = scorpionInput$gex,
                           ppiNet = scorpionInput$ppi)
meltedScorpionOutput <-melt(as.matrix(scorpionOutput$regNet))
colnames(meltedScorpionOutput) <- c("tf", "gene", "score")
write.csv(meltedScorpionOutput, paste0(args[2], ".csv"))
