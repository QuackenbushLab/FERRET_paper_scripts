# Dependency on Slingshot for pseudotime inference
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("monocle3")){
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'HDF5Array',
                         'terra', 'ggrastr'))
  devtools::install_github('cole-trapnell-lab/monocle3')
}
pseudotime_result_dir <- NULL
if(!file.exists(paste0(pseudotime_result_dir, "/final_assignment/"))){
  dir.create(paste0(pseudotime_result_dir, "/final_assignment/"))
}
for(file in file_subset){
    monocle_expressions_graph <- readRDS(paste0(pseudotime_result_dir, "/ptime_", file, ".RDS"))
    
    # Assign pseudotime trajectory.
    monocle_expressions_pseudotime <- monocle3::order_cells(monocle_expressions_graph)
    plot_cells(monocle_expressions_pseudotime, color_cells_by = "pseudotime")
    ggsave(paste0(pseudotime_result_dir, "/final_assignment/plot_", file, ".pdf"), plot = last_plot()) 
    pseudotime <- monocle3::pseudotime(monocle_expressions_pseudotime)
    write.csv(pseudotime, paste0(pseudotime_result_dir, "/final_assignment/pseudotime_", file)) 
}