# Set a random seed.
set.seed(1)

# Split on only the data that have a clear pseudotime trajectory.
file_subset <- c("1188_Ru1170gFreeze-sort_IGO_10034_8_dense.csv_Endothelial cells_logcounts.csv",
                 "1188_Ru1170gFreeze-sort_IGO_10034_8_dense.csv_Macrophages_logcounts.csv",
                 "1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv_Endothelial cells_logcounts.csv",
                 "1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv_Fibroblasts_logcounts.csv",
                 "1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv_Macrophages_logcounts.csv",
                 "1189_Ru1170gsort-freeze_IGO_10034_10_dense.csv_NK cells_logcounts.csv",
                 "1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv_B cells_logcounts.csv",
                 "1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv_Fibroblasts_logcounts.csv",
                 "1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv_Macrophages_logcounts.csv",
                 "1202_Ru1134_Freeze_thaw_IGO_10034_12_dense.csv_T cells_logcounts.csv",
                 "1663_Ru1271_IGO_10414_11_dense.csv_Airway goblet cells_logcounts.csv",
                 "1663_Ru1271_IGO_10414_11_dense.csv_Macrophages_logcounts.csv",
                 "RU661_TUMOR_dense.csv_Macrophages_logcounts.csv",
                 "RU661_TUMOR_dense.csv_T cells_logcounts.csv",
                 "RU682_NORMAL_dense.csv_Macrophages_logcounts.csv",
                 "RU682_NORMAL_dense.csv_NK cells_logcounts.csv",
                 "RU682_NORMAL_dense.csv_T cells_logcounts.csv",
                 "RU682_TUMOR_dense.csv_Endothelial cells_logcounts.csv",
                 "RU682_TUMOR_dense.csv_Fibroblasts_logcounts.csv",
                 "RU682_TUMOR_dense.csv_Macrophages_logcounts.csv",
                 "RU682_TUMOR_dense.csv_T cells_logcounts.csv",
                 "RU684_NORMAL_dense.csv_Macrophages_logcounts.csv",
                 "RU684_NORMAL_dense.csv_T cells_logcounts.csv",
                 "Ru1135_dense.csv_Fibroblasts_logcounts.csv",
                 "Ru1135_dense.csv_Macrophages_logcounts.csv",
                 "Ru1137_dense.csv_Endothelial cells_logcounts.csv",
                 "Ru1137_dense.csv_Macrophages_logcounts.csv")

# Directories
split_result_dir <- NULL
expression_dir <- NULL
pseudotime_dir <- NULL
pseudotime_result_dir <- NULL
if(!file.exists(split_result_dir)){
  dir.create(split_result_dir)
  dir.create(pseudotime_result_dir)
}

# For each file, create a 5-fold cross-validation. For each, 
# Split the pseudotime trajectory and the expression data.
# Also, copy the original pseudotime trajectory and expression data.
for(file in file_subset){
  
  # Read the expression and pseudotime files.
  expression <-  read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                                           header = TRUE, row.names = 1)
  pseudotime <- read.csv(paste0(pseudotime_dir, "/pseudotime_", 
                                  file), header = TRUE, row.names = 1)
  rownames(pseudotime) <- colnames(expression)

  write.table(expression, paste0(split_result_dir, "/", file, "_original.tsv"), sep = "\t",
              quote = FALSE, row.names = TRUE, col.names = NA)
  write.csv(pseudotime, paste0(pseudotime_result_dir, "/pseudotime_", file, "_original.csv"))

  # Do 1 random 5-fold split.
  folds <- 1
  # Randomly permute the cells.
  perm <- sample(colnames(expression), ncol(expression), replace = FALSE)
  expression_perm <- expression[,perm]
  pseudotime_perm <- as.data.frame(pseudotime[perm, ])
  rownames(pseudotime_perm) <- perm
  colnames(pseudotime_perm) <- "x"

  # Perform the split.
  for(j in 1:folds){

    # Subset everything not in the fold.
    fold_size <- ceiling(ncol(expression_perm) / 5)
    start_pos <- (j-1) * fold_size + 1
    end_pos <- j * fold_size
    print(paste(i, j, start_pos, end_pos))
    to_keep <- setdiff(colnames(expression_perm), colnames(expression_perm)[start_pos:end_pos])

    # Save everything not in the fold.
    expressions_to_keep <- expression_perm[,to_keep]
    pseudotime_to_keep <- as.data.frame(pseudotime_perm[to_keep,])
    rownames(pseudotime_to_keep) <- to_keep
    colnames(pseudotime_to_keep) <- "x"
    write.table(expressions_to_keep, paste0(split_result_dir, "/", file, "_", i, "_", j, ".tsv"), sep = "\t",
                quote = FALSE, row.names = TRUE, col.names = NA)
    write.csv(pseudotime_to_keep, paste0(pseudotime_result_dir, "/pseudotime_", file, "_", i, "_", j, ".csv"))
  }
}