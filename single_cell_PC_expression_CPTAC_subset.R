if(!require(irlba)){
  install.packages(irlba)
}
library(irlba)

# Split on only the data that have a clear pseudotime trajectory.
set.seed(1)
sampNames <- c("C3L-00359-01", "C3N-00662-03", "C3N-01814-01", "C3N-01904-02",
               "C3N-02181-02", "C3N-02784-01")
files <- list(c("C3L-00359-01_Macrophages", "C3L-00359-01_Podocytes"),
           c("C3N-00662-03_Neurons", "C3N-00662-03_Oligodendrocytes"),
           c("C3N-01814-01_Interneurons", "C3N-01814-01_Macrophages",
           "C3N-01814-01_Neurons", "C3N-01814-01_Oligodendrocytes"),
           c("C3N-01904-02_Macrophages", "C3N-01904-02_Podocytes"),
           c("C3N-02181-02_Interneurons", "C3N-02181-02_Macrophages",
           "C3N-02181-02_Neurons", "C3N-02181-02_Oligodendrocytes"),
           c("C3N-02784-01_Macrophages", "C3N-02784-01_Oligodendrocytes"))

# Directories
pca_result_dir <- NULL
plot_result_dir <- NULL
expression_dir <- NULL
split_dir <- NULL

# For each sample, get the whole "counts" file and do a PCA, then find the elbow
# point.
for(sample in sampNames){
  
  # Read the original count file.
  expression <-  t(read.csv(paste(expression_dir, paste0(sample, "_counts.csv"), sep = "/"),
                                           header = TRUE, row.names = 1))
  
  # Do PCA.
  pca <- irlba::prcomp_irlba(expression, n = 100)
  saveRDS(pca, paste0(plot_result_dir, "/", sample, "_pca.RDS"))
  
}

# For all files, project onto PC's.
for(i in 1:length(sampNames)){
  
  # Read PCA.
  pca <- readRDS(paste0(plot_result_dir, "/", sampNames[i], "_pca.RDS"))

  # Project all files.
  for(cellType in files[[i]]){
    # Read the expression file.
    expression <- t(read.table(paste0(split_dir, "/", cellType, "_original.tsv"), sep = "\t",
                               header = TRUE, row.names = 1))
    
    # Project.
    expPCA <- expression %*% pca$rotation
    
    # Save.
    write.table(expPCA, paste0(pca_result_dir, "/", cellType, "_original_pca.tsv"))
    
    for(fold in 1:10){
      for(iteration in 1:5){
        
        # Read the expression file.
        expression <- t(read.table(paste0(split_dir, "/", cellType, "_", iteration, "_", fold, ".tsv"), sep = "\t",
                    header = TRUE, row.names = 1))
        
        # Project.
        expPCA <- expression %*% pca$rotation
        
        # Save.
        write.table(expPCA, paste0(pca_result_dir, "/", cellType, "_", iteration, "_", fold, "_pca.tsv"))
      }
    }
  }
}