# Directories
split_dir <- NULL
transpose_split_dir <- NULL

# For all files, project onto PC's.
for(i in list.files(split_dir)){
    expPCA <- read.table(paste0(split_dir, "/", i), sep = " ")
    write.table(t(expPCA), paste0(transpose_split_dir, "/", i), sep = " ", col.names = NA, row.names = TRUE)
}
