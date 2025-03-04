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

# This function obtains the position of each cell in the PDF of
# each gene of interest, then averages them over all genes.
obtain_avg_density <- function(exp){
  # Get CDF.
  densities <- cbind(lapply(rownames(exp), function(gene){
    density <- density(as.matrix(exp[gene,]))
    pdf <- density$y / sum(density$y)
    cdf <- cumsum(pdf)
    names(cdf) <- density$x
    return(cdf)
  }))
  # For each cell, find the average percentile across all genes using CDF.
  avg_density_per_cell <- unlist(lapply(1:ncol(exp), function(i){
    all_d <- unlist(lapply(1:length(densities), function(d){
      closest_index <- which.min(abs(as.numeric(names(densities[[d]])) - exp[d,i]))
      density <- densities[[d]][closest_index]
      return(density)
    }))
    return(mean(all_d))
  }))
  return(avg_density_per_cell)
}

# This method computes the UMAP, clusters, and graph and plots them with respect
# to the genes of interest.
find_graph <- function(differentiationTraj, expressions, cell_size){
  # Read in trajectory.
  metadata <- data.frame(estTrajectory = differentiationTraj)
  rownames(metadata) <- colnames(expressions)
  colnames(metadata) <- "estTrajectory"
  expressionsMat <- as.matrix(expressions)
  
  # Create data set.
  monocle_expressions <- monocle3::new_cell_data_set(expressionsMat, cell_metadata = metadata)
  
  # Wrangle trajectory so that we can use it for plotting.
  names(monocle_expressions@colData$estTrajectory) <- colnames(expressions)

  # Reduce dimensionality using UMAP. Running LSI is a precursor to this.
  monocle_expressions_pca <- monocle3::preprocess_cds(monocle_expressions, method = "LSI")
  monocle_expressions_umap <- monocle3::reduce_dimension(monocle_expressions_pca, reduction_method = "UMAP",
                                                         preprocess_method = "LSI")
  # Clustering the cells using the UMAP input.
  monocle_expressions_cluster <- monocle3::cluster_cells(monocle_expressions_umap, reduction_method = "UMAP")
  print("finished clustering")
  monocle_expressions_graph <- monocle3::learn_graph(monocle_expressions_cluster)
  saveRDS(monocle_expressions_graph, paste0(pseudotime_result_dir, "/ptime_", file, ".RDS"))
  plot_cells(monocle_expressions_graph, color_cells_by = "estTrajectory", cell_size = cell_size)
  ggsave(paste0(pseudotime_result_dir, "/plotTraj_", file, ".pdf"), plot = last_plot())
}

# For macrophages, combine the markers of M1 and M2 macrophages to designate
# a trajectory from M0 -> M1 -> M2.
# These markers are taken from Boutilier et al.
# We mark the time on a gradation from 0 (M0) to 1 (M1) to 2 (M2).
M0 <- unname(mapIds(org.Hs.eg.db,
             keys=c("CSF1R", "CD14", "CD68", "CD11B"),
             column="ENSEMBL",
             keytype="SYMBOL",
             multiVals="first"))
M0 <- M0[which(!is.na(M0))]
M1 <- unname(mapIds(org.Hs.eg.db,
                    keys=c("CD86", "MARCO", "CXCL9", "CXCL10", "CXCL11", "NOS2", "SOCS1", "CD64"),
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first"))
M1 <- M1[which(!is.na(M1))]
M2 <- unname(mapIds(org.Hs.eg.db,
                    keys=c("TGM2", "CD23", "ARG1", "CCL22", "CD163", "CD206"),
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first"))
M2 <- M2[which(!is.na(M2))]
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Macrophages_counts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.table(paste(expression_dir, file, sep = "/"), sep = "\t",
                              header = TRUE, row.names = 1)
    
    # Obtain a density distribution for each of the types.
    m0 <- intersect(M0, rownames(expressions))
    m1 <- intersect(M1, rownames(expressions))
    m2 <- intersect(M2, rownames(expressions))
    M0_distrib <- obtain_avg_density(expressions[m0,])
    M1_distrib <- obtain_avg_density(expressions[m1,])
    M2_distrib <- obtain_avg_density(expressions[m2,])
    distrib <- data.frame(m0 = M0_distrib, m1 = M1_distrib, m2 = M2_distrib)
    
    # Find the max across all distributions.
    which_max <- apply(distrib, 1, which.max)
    max_2_3 <- apply(distrib[,2:3], 1, max)
    max_1_3 <- apply(distrib[,c(1,3)], 1, max)
    max_1_2 <- apply(distrib[,1:2], 1, max)
    max <- apply(distrib, 1, max)
    
    # Quantify the extent to which a cell is of one type vs. another.
    metric <- rep(NA, nrow(distrib))
    cutoff <- 0.1
    which_m0_meets_cutoff <- intersect(which(which_max == 1), which(max - max_2_3 > cutoff))
    which_m1_meets_cutoff <- intersect(which(which_max == 2), which(max - max_1_3 > cutoff))
    which_m2_meets_cutoff <- intersect(which(which_max == 3), which(max - max_1_2 > cutoff))
    metric[which_m0_meets_cutoff] <- 0.5 - ((max[which_m0_meets_cutoff] - M1_distrib[which_m0_meets_cutoff]) / 2)
    metric[which_m1_meets_cutoff] <- 1 - ((2 * max[which_m1_meets_cutoff] - M2_distrib[which_m1_meets_cutoff]
                                          - M0_distrib[which_m1_meets_cutoff]) / 4)
    metric[which_m2_meets_cutoff] <- 1.5 + ((max[which_m2_meets_cutoff] - M1_distrib[which_m2_meets_cutoff]) / 2)
    pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
    hist(metric)
    dev.off()
    print(dim(expressions))
    print(length(metric))
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    find_graph(metric, expressions, 1) 
  }
}
names(time_list_macrophages) <- list.files(expression_dir)
time_list_macrophages <- which(!is.null(time_list_macrophages))

# For podocytes, LOX is prognostic of worse survival (Chen et al, Anazco et al).
# EGFR is prognostic of worse survival (Moch, Weber, Ciardello, Minner)
# DACH1 is essential for podocyte function and is lost during dedifferentiation (Endlich et al).
# GLEPP1 (PTPRO) decreases during dedifferentiation (Sharif et al)
podo_worse <- unname(mapIds(org.Hs.eg.db,
                            #keys=c("LOX", "EGFR"), # Not related to UMAP
                            #keys="EGFR", # Not related to UMAP
                            #keys="WT1", # Doesn't overlap with expression data
                            #keys="NOTCH3", # Not enough variation
                            #keys="PTPRO", # Not enough variation
                            keys=c("DACH1", "PTPRO"),
                            column="ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first"))

# For each file, plot the distribution of the key genes.
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Podocytes_counts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.table(paste(expression_dir, file, sep = "/"), sep = "\t",
                              header = TRUE, row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(podo_worse, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
    hist(metric)
    dev.off()
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    find_graph(metric, expressions, 1) 
  }
}

# In neurons, NOTCH1 and HES1 are overexpressed with increased proliferation
# of tumor cells. We consider a trajectory from less proliferative to
# more proliferative. (Cenciarelli et al)
npc <- unname(mapIds(org.Hs.eg.db,
                            #keys="SOX2", # No relationship to UMAP
                            #keys = "HES1", # No relationship to UMAP
                            keys = c("HES1", "NOTCH1"),
                            column="ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first"))

# For each file, plot the distribution of the key genes.
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Neurons_counts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.table(paste(expression_dir, file, sep = "/"), sep = "\t",
                              header = TRUE, row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(npc, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
    hist(metric)
    dev.off()
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    find_graph(metric, expressions, 0.75) 
  }
}
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Interneurons_counts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.table(paste(expression_dir, file, sep = "/"), sep = "\t",
                              header = TRUE, row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(npc, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
    hist(metric)
    dev.off()
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    find_graph(metric, expressions, 0.75) 
  }
}

# In oligodendrocytes, myelin components are under-expressed in gliomas
# (at least for IDH-mut), suggesting stalled differentiation (Wei et al).
npc <- unname(mapIds(org.Hs.eg.db,
                     keys = c("MBP", "MOG"),
                     column="ENSEMBL",
                     keytype="SYMBOL",
                     multiVals="first"))

# For each file, plot the distribution of the key genes.
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Oligodendrocytes_counts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.table(paste(expression_dir, file, sep = "/"), sep = "\t",
                              header = TRUE, row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(npc, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
    hist(metric)
    dev.off()
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    find_graph(metric, expressions, 0.75) 
  }
}