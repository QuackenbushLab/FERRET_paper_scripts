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
# These packages need to be installed so that as_cholmod_sparse works in the call to preprocess_cds.
# Note that we need to run brew install gcc prior to installing Matrix. Then,
# create a new file called ~/.R/Makevars and add the following lines:
# FC = /opt/homebrew/Cellar/gcc/13.2.0/bin/gfortran
# F77 = /opt/homebrew/Cellar/gcc/13.2.0/bin/gfortran
# FLIBS = -L//opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/13
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
library("ggplot2")
library("org.Hs.eg.db")

# Directories
expression_dir <- NULL
pseudotime_result_dir <- NULL

if (!file.exists(pseudotime_result_dir)){
  dir.create(file.path(pseudotime_result_dir))
}
 
# For each file, find a pseudotime trajectory.
if(!file.exists(paste0(pseudotime_result_dir, "/final_assignment/"))){
  dir.create(paste0(pseudotime_result_dir, "/final_assignment/"))
}
for(i in 1:length(list.files(expression_dir))){
  file <- list.files(expression_dir)[i]
  if(length(grep("s_counts", file)) > 0){
  #if(length(grep("Podocytes_counts", file)) > 0){
    if(!file.exists(paste0(pseudotime_result_dir, "/final_assignment/pseudotime_", file))){
      if(file.exists(paste0(pseudotime_result_dir, "/ptime_", file, ".RDS"))){
        monocle_expressions_graph <- readRDS(paste0("/Users/tae771/Documents/postdoc/single_cell_TCGA/pseudotime/ptime_", file, ".RDS"))
        
        # Assign pseudotime trajectory.
        monocle_expressions_pseudotime <- monocle3::order_cells(monocle_expressions_graph)
        plot_cells(monocle_expressions_pseudotime, color_cells_by = "pseudotime")
        ggsave(paste0("/Users/tae771/Documents/postdoc/single_cell_TCGA/pseudotime/final_assignment/plot_", file, ".pdf"), plot = last_plot()) 
        pseudotime <- monocle3::pseudotime(monocle_expressions_pseudotime)
        write.csv(pseudotime, paste0("/Users/tae771/Documents/postdoc/single_cell_TCGA/pseudotime/final_assignment/pseudotime_", file)) 
      }
    }
  }
}

# Re-do several of the pseudotime trajectories.
pseudotime_redo <- c(#"C3N-00662-03_Oligodendrocytes_counts.tsv",
                     #"C3N-01814-01_Interneurons_counts.tsv",
                     #"C3N-01904-02_Podocytes_counts.tsv",
                     "C3N-02181-02_Interneurons_counts.tsv",
                     "C3N-02181-02_Oligodendrocytes_counts.tsv",
                     "C3N-02784-01_Macrophages_counts.tsv")
for(file in pseudotime_redo){
  monocle_expressions_graph <- readRDS(paste0(pseudotime_result_dir, "/ptime_", file, ".RDS"))
  
  # Assign pseudotime trajectory.
  monocle_expressions_pseudotime <- monocle3::order_cells(monocle_expressions_graph)
  plot_cells(monocle_expressions_pseudotime, color_cells_by = "pseudotime")
  ggsave(paste0(pseudotime_result_dir, "/final_assignment/plot_", file, ".pdf"), plot = last_plot()) 
  pseudotime <- monocle3::pseudotime(monocle_expressions_pseudotime)
  write.csv(pseudotime, paste0(pseudotime_result_dir, "/final_assignment/pseudotime_", file))
}

# For neurons in C3N-00662-03_Neurons_counts.tsv, make the range of values
# from UMAP X = c(-8,0) Y = c(-8,-3) range from 0 to 18.
monocle_expressions_graph <- readRDS(paste0(pseudotime_result_dir, "/ptime_C3N-00662-03_Neurons_counts.tsv.RDS"))
umap <- monocle_expressions_graph@reduce_dim_aux@listData$UMAP@listData$model@listData$umap_model$embedding
positions_lower_left_quadrant <- rownames(umap)[Reduce(intersect, list(which(umap[,1] > -8),
                                                                       which(umap[,1] < 0),
                                                                       which(umap[,2] > -8),
                                                                       which(umap[,2] < -3)))]
monocle_expressions_pseudotime <- read.csv(paste0(pseudotime_result_dir, 
                                                    "/final_assignment/pseudotime_C3N-00662-03_Neurons_counts.tsv"),
                                           row.names = 1)
range_lower_left_quadrant <- range(monocle_expressions_pseudotime[positions_lower_left_quadrant,"x"])
factor <- 18 / range_lower_left_quadrant[2]
scaled_lower_left_quadrant <- monocle_expressions_pseudotime[positions_lower_left_quadrant,"x"] * factor
monocle_expressions_graph@principal_graph_aux@listData$UMAP$pseudotime <- monocle_expressions_pseudotime$x
names(monocle_expressions_graph@principal_graph_aux@listData$UMAP$pseudotime) <- rownames(monocle_expressions_pseudotime)
monocle_expressions_graph@principal_graph_aux@listData$UMAP$pseudotime[positions_lower_left_quadrant] <- monocle_expressions_pseudotime[positions_lower_left_quadrant,"x"] * factor
plot_cells(monocle_expressions_graph, color_cells_by = "pseudotime")
ggsave(paste0(pseudotime_result_dir, "/final_assignment/plot_C3N-00662-03_Neurons_counts_manual.pdf"), plot = last_plot()) 
pseudotime <- monocle3::pseudotime(monocle_expressions_graph)
write.csv(pseudotime, paste0(pseudotime_result_dir, "/final_assignment/pseudotime_", file))