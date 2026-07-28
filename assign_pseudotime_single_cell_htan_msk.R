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
if(!require("Matrix")){
  install.packages("Matrix", type = "source")
}
if(!require("irlba")){
  install.packages("irlba", type = "source")
}
library("ggplot2")
library("org.Hs.eg.db")

# Directories
dataDir <- "/Users/tae771/Documents/postdoc/single_cell_HTAN_MSK/"
expression_dir <- paste0(dataDir, "/cellTypeSpecificMatrices")
pseudotime_result_dir <- paste0(dataDir, "/pseudotime")

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
M0 <- c("CSF1R", "CD14", "CD68", "CD11B")
M1 <- c("CD86", "MARCO", "CXCL9", "CXCL10", "CXCL11", "NOS2", "SOCS1", "CD64")
M2 <- c("TGM2", "CD23", "ARG1", "CCL22", "CD163", "CD206")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Macrophages_logcounts", file)) > 0){

    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), row.names = 1)
    
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
    tryCatch({
      find_graph(metric, expressions, 1) 
    }, error = function(cond){
      print(paste(file, "was too small!"))
    })
    
  }
}

# In B cell progenitors, E2A and HEB induce B-lineage specific gene expression.
# These genes are from Miyazaki et al.
B_differentiation <- c("TCF3", "FOXO1", "EBF1", "PAX5")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("B cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                              row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(B_differentiation, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# Genes upregulated in the T cell regulatory network are also documented in
# Miyazaki et al
T_differentiation <- c("TCF3", "NOTCH1", "HES1", "TCF7", "TLE3", "TLE6", "BCL11B",
                       "GATA3")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("T cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(T_differentiation, rownames(expressions))
    str(pw)
    #pw_distrib <- obtain_avg_density(expressions[pw,])
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# Clara cell differentiation is documented in Ustiyan et al.
Clara_differentiation <- "FOXM1"
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Clara cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(Clara_differentiation, rownames(expressions))
    str(pw)
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# Endothelial cell differentiation is documented in Kanki et al.
Endothelial_differentiation <- c("GATA2", "FLI1", "SOX7", "SOX18")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Endothelial cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(Endothelial_differentiation, rownames(expressions))
    str(pw)
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# For fibroblasts, we look at the process of fibrogenesis (which occurs in
# inflammatory response). This is documented in Darby and Hewitson.
Fibrogenesis <- c("TGFB1", "FGF2", "ANGPT1", "ANGPT2", "PDGFA", "PDGFB",
                  "PDGFC", "PDGFD", "VEGFA", "VEGFB")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Fibroblasts_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(Fibrogenesis, rownames(expressions))
    str(pw)
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# Pericytes are precursors to myofibroblasts, so we can also

# Genes involved in mast cell development and survival are discussed in Jayapal et al.
Mast_development <- c("KITLG", "IL3")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Mast cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(Mast_development, rownames(expressions))
    str(pw)
    metric <- as.numeric(colMeans(expressions[pw,]))
    tryCatch({
      pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
      hist(metric)
      dev.off()
    }, error = function(cond){
      print(cond)
    })
    
    write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                col.names = FALSE, row.names = FALSE, sep = "\n")
    
    tryCatch({
      find_graph(metric, expressions, 1)
    }, error = function(cond){
      print(paste(file, "is too small!"))
    })
  }
}

# Genes involved in ROS formation and NET production in neutrophils (both part of
# immune response) are described in Burn et al.
ros_net <- c("G6PD", "TKT")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Neutrophils_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(ros_net, rownames(expressions))
    
    tryCatch({
        metric <- as.numeric(colMeans(expressions[pw,]))
      tryCatch({
        pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
        hist(metric)
        dev.off()
      }, error = function(cond){
        print(cond)
      })
      
      write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                  col.names = FALSE, row.names = FALSE, sep = "\n")
      
      tryCatch({
        find_graph(metric, expressions, 1)
      }, error = function(cond){
        print(paste(file, "is too small!"))
      })
    }, error = function(cond){
      print("only one cell!")
      str(expressions[pw,])
    })
  }
}

# In NK cells, FASL and TRAIL are transmembrane proteins that mediate cytotoxicity.
# (Prager and Watzl)
nk_cytotoxicity <- c("FASLG", "TNFSF10")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("NK cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(nk_cytotoxicity, rownames(expressions))
    
    tryCatch({
      metric <- as.numeric(colMeans(expressions[pw,]))
      tryCatch({
        pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
        hist(metric)
        dev.off()
      }, error = function(cond){
        print(cond)
      })
      
      write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                  col.names = FALSE, row.names = FALSE, sep = "\n")
      
      tryCatch({
        find_graph(metric, expressions, 1)
      }, error = function(cond){
        print(paste(file, "is too small!"))
      })
    }, error = function(cond){
      print("only one cell!")
      str(expressions[pw,])
    })
  }
}

# We skip Pulmonary Alveolar Type II cells as their markers are not well-studied.

# The ependymal cells are likely cancerous epithelial cells, not ependymal cells.
# We leave these out for now. However, they have been previously found in metaplastic
# carcinoma. Bard cites the following. However, I cannot find these papers.
# 1. "Lung adenocarcinoma with extensive squamous and glandular differentiation 
# featuring a prominent component of ependymal-like cells: A case report and literature 
# review" by Tsuchida et al. (2014)
# # 
# # 2. "Spindle and ependymal-like carcinoma of the lung: A case report and literature 
# review" by Shen et al. (2017) 
# # 
# # 3. "Ependymal-like differentiation in pulmonary adenocarcinoma: A morphologic 
# and immunohistochemical analysis of four cases" by Pelosi et al. (2008)
# # 
# # 4. "Atypical pulmonary carcinomas with ependymal-like differentiation: Report of 
# two cases and review of the literature" by Nakamura et al. (2005) 
# # 
# # 5. "Ependymal-like differentiation in a pulmonary adenocarcinoma: Clinicopathologic 
# and ultrastructural observations" by Koga et al. (2002)

# Could not find a clear trajectory for pericytes.

# Several markers of goblet cells have been increased by exposing HAECs to IL-17.
# (Rahmawati et al)
goblet_differentiation <- c("MUC5AC", "MUC5B", "SPDEF")
for(file in list.files(expression_dir)){
  metric <- NULL
  if(length(grep("Airway goblet cells_logcounts", file)) > 0){
    print(file)
    # Read the file.
    expressions <- read.csv(paste(expression_dir, file, sep = "/"), sep = ",",
                            row.names = 1)
    
    # Obtain a density distribution for each of the types.
    pw <- intersect(goblet_differentiation, rownames(expressions))
    
    tryCatch({
      metric <- as.numeric(colMeans(expressions[pw,]))
      tryCatch({
        pdf(paste0(pseudotime_result_dir, "/hist_", file, ".pdf"))
        hist(metric)
        dev.off()
      }, error = function(cond){
        print(cond)
      })
      
      write.table(metric, paste0(pseudotime_result_dir, "/ptime_", file, ".csv"),
                  col.names = FALSE, row.names = FALSE, sep = "\n")
      
      tryCatch({
        find_graph(metric, expressions, 1)
      }, error = function(cond){
        print(paste(file, "is too small!"))
      })
    }, error = function(cond){
      print("only one cell!")
      str(expressions[pw,])
    })
  }
}

 
# For each file, find a pseudotime trajectory. Subset the files that actually
# showed a meaningful trajectory.
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
file_subset <- c("RU682_NORMAL_dense.csv_Macrophages_logcounts.csv")
if(!file.exists(paste0(pseudotime_result_dir, "/final_assignment/"))){
  dir.create(paste0(pseudotime_result_dir, "/final_assignment/"))
}