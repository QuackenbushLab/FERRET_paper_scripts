# Split on only the data that have a clear pseudotime trajectory.
files <- c("C3L-00359-01_Macrophages", "C3L-00359-01_Podocytes",
           "C3N-00662-03_Neurons", "C3N-00662-03_Oligodendrocytes",
           "C3N-01814-01_Interneurons", "C3N-01814-01_Macrophages",
           "C3N-01814-01_Neurons", "C3N-01814-01_Oligodendrocytes",
           "C3N-01904-02_Macrophages", "C3N-01904-02_Podocytes",
           "C3N-02181-02_Interneurons", "C3N-02181-02_Macrophages",
           "C3N-02181-02_Neurons", "C3N-02181-02_Oligodendrocytes",
           "C3N-02784-01_Macrophages", "C3N-02784-01_Oligodendrocytes")

# Directories
split_result_dir <- NULL
expression_dir <- NULL
split_result_dir_pseudotime <- NULL
if(!file.exists(split_result_dir_pseudotime)){
  dir.create(split_result_dir_pseudotime)
}

# For each file, take a random subset of 63 cells and create 5 10-fold cross-validations. For each, 
# Split the pseudotime trajectory and the expression data.
# Also, copy the original pseudotime trajectory and expression data.
for(file in files){
  
  # Read the expression and pseudotime files.
  expression <-  read.table(paste0(split_result_dir, "/", file, "_original.tsv"),
                            sep = "\t", header = TRUE, row.names = 1)
  expression_unfiltered <- read.table(paste(expression_dir, paste0(file, "_counts.tsv"), sep = "/"), sep = "\t",
                                                     header = TRUE, row.names = 1)
  pseudotime <- read.csv(paste0(pseudotime_result_dir, "/pseudotime_", 
                                  file, "_counts.tsv"), header = TRUE, row.names = 1)
  rownames(pseudotime) <- colnames(expression_unfiltered)

  # Subset.
  pseudotime_subset <- as.data.frame(pseudotime[colnames(expression),])
  rownames(pseudotime_subset) <- colnames(expression)
  colnames(pseudotime_subset) <- "x"
  
  write.csv(pseudotime_subset, paste0(split_result_dir_pseudotime, "/pseudotime_", file, "_original.csv"))

  # Do 5 random 10-fold splits.
  split_count <- 5
  folds <- 10
  for(i in 1:split_count){

    # Perform the split.
    for(j in 1:folds){
      # Read expression data.
      expression_to_keep <- read.table(paste0(split_result_dir, "/", file, "_", i, "_", j, ".tsv"),
                                       sep = "\t", header = TRUE, row.names = 1)
      print(dim(expression_to_keep))
      pseudotime_to_keep <- as.data.frame(pseudotime_subset[colnames(expression_to_keep),])
      rownames(pseudotime_to_keep) <- colnames(expression_to_keep)
      colnames(pseudotime_to_keep) <- "x"
      str(pseudotime_to_keep)
      str(rownames(pseudotime_to_keep))
      write.csv(pseudotime_to_keep, paste0(split_result_dir_pseudotime, "/pseudotime_", file, "_", i, "_", j, ".csv"))
    }
  }
}