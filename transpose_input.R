# Directories
split_dir <- "/home/ubuntu/CPTAC_single_cell_splits_subset_pca"
transpose_split_dir <- "/home/ubuntu/CPTAC_single_cell_splits_subset_pca_t"

# For all files, project onto PC's.
for(i in list.files(split_dir)){
    expPCA <- read.table(paste0(split_dir, "/", i), sep = " ")
    write.table(t(expPCA), paste0(transpose_split_dir, "/", i), sep = " ", col.names = NA, row.names = TRUE)
}
