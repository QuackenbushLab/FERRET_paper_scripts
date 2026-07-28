# Need reshape2 to create the adjacency list.
library(reshape2)

# Read in all cutoffs for each PIDC result file, then format them for FERRET.
sourceFileList <- c(rep("C3L-00359-01", 2), 
                    rep("C3N-00662-03", 2), 
                    rep("C3N-01814-01", 4),
                    rep("C3N-01904-02", 2),
                    rep("C3N-02181-02", 4),
                    rep("C3N-02784-01", 2))
ingroupName <- c("Macrophages", "Podocytes", 
                 "Oligodendrocytes", "Neurons", 
                 "Interneurons", "Macrophages", "Neurons","Oligodendrocytes",
                 "Macrophages", "Podocytes",
                 "Interneurons", "Macrophages", "Neurons", "Oligodendrocytes",
                 "Macrophages", "Oligodendrocytes")
cutoffs <- rev(seq(0.1, 0.9, 0.1))

# Create new directory.
originalDir <- NULL
consolidatedDir <- NULL
if(!file.exists(consolidatedDir)){
  dir.create(consolidatedDir)
}

# Loop through files.
for(i in 1:length(sourceFileList)){
  
  # Initialize the first adjacency matrix.
  adjacency <- matrix(rep(0, 100 * 100), nrow = 100)
  rownames(adjacency) <- paste0("PC", 1:100)
  colnames(adjacency) <- rownames(adjacency)
  
  # Loop through the randomizations.
  for(j in 1:5){
    # Loop through the folds.
    for(k in 1:10){
      # Loop through the cutoffs backwards.
      for(cutoff in cutoffs){
        
        # Read the file.
        adjacencyCutoff <- read.csv(file = paste0(originalDir, "/", cutoff, "_", sourceFileList[i],
                                                  "_", ingroupName[i], "_", j, "_", k, "_pca.tsv"))
        adjacency[which(adjacencyCutoff == "true")] <- 1 - cutoff
      }
      
      # Create adjacency list from adjacency matrix.
      adjacencyList <- melt(adjacency)
      names(adjacencyList) <- c("source", "target", "score")
      adjacencyListNonzero <- adjacencyList[which(adjacencyList$score > 0),]
      write.csv(adjacencyListNonzero, file = paste0(consolidatedDir, "/", sourceFileList[i],
                              "_", ingroupName[i], "_", j, "_", k, "_pca.tsv"))
    }
  }
  
  # Next, do originals.
  for(cutoff in cutoffs){
    
    # Read the file.
    adjacencyCutoff <- read.csv(file = paste0(originalDir, "/", cutoff, "_", sourceFileList[i],
                                              "_", ingroupName[i], "_original_pca.tsv"))
    adjacency[which(adjacencyCutoff == "true")] <- 1 - cutoff
  }
  
  # Create adjacency list from adjacency matrix.
  adjacencyList <- melt(adjacency)
  names(adjacencyList) <- c("source", "target", "score")
  adjacencyListNonzero <- adjacencyList[which(adjacencyList$score > 0),]
  write.csv(adjacencyListNonzero, file = paste0(consolidatedDir, "/", sourceFileList[i],
                                                "_", ingroupName[i], "_original_pca.tsv"))
}
