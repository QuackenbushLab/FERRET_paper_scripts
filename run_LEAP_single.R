# Reference: LEAP vignette (http://cran.nexr.com/web/packages/LEAP/vignettes/LEAP_Vignette.pdf)
library("LEAP")
library("reshape2")

# Load arguments.
args = commandArgs(trailingOnly = TRUE)

# Load files.
expression <- read.table(args[1], sep = " ", header = TRUE, row.names = 1)
pseudotime <- read.csv(args[2], row.names = 1)

# Rename columns of expression matrix and sort by pseudotime.
ptime_names <- make.names(pseudotime[colnames(expression), "x"])
expression_ptime <- expression
colnames(expression_ptime) <- ptime_names
expression_ptime_sorted <- expression_ptime[,order(pseudotime$x)]

# Compute MAC. We want to retain as much information as possible,
# so set the MAC cutoff to 0. We set max_lag_prop to 1/3 as suggested
# by the vignette.
MAC_results <- MAC_counter(data = expression_ptime_sorted, max_lag_prop = 1/3,
			   file_name = FALSE, MAC_cutoff = 0)

# Compute the significance of results using permutation. Use 100 permutations
# as suggested by the vignette. We set FDR_cutoffs to 101 as done in the vignette.
MAC_perm_results <- MAC_perm(data = expression_ptime_sorted, MACs_observ = MAC_results)

# Substitute FDR scores for correlations in the MAC_results.
MAC_cor_rounded <- unlist(lapply(unname(MAC_results[,"Correlation"]), function(cor){return(abs(round(cor, 2)))}))
scored_adj_list <- data.frame(source = paste0("PC",MAC_results[,"Row gene index"]),
			      target = paste0("PC",MAC_results[,"Column gene index"]),
			      score = unlist(lapply(MAC_cor_rounded, function(cor){return(MAC_perm_results[which.min(abs(MAC_perm_results[,"cors"] - cor)), "fdr"])})))
write.csv(scored_adj_list, args[3])
