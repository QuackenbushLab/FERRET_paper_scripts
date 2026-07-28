if(!require("R.matlab")){
	install.packages("R.matlab")
}
library("R.matlab")

expressionPath = NULL
pseudotimePath = NULL
outputPath = NULL
dir.create(outputPath)

# Loop through each file and process it
for(file in list.files(expressionPath)){
	# Read expression data.
	expression = as.matrix(read.table(paste0(expressionPath, "/", file), header = TRUE, sep = " ", row.names = 1))

   	# Read the pseudotime file.
   	pseudotime = t(as.matrix(read.csv(paste0(pseudotimePath, "/pseudotime_", strsplit(file, "_pca.tsv")[[1]][1], ".csv"), row.names = 1)))

	# Save the data.
	str(file)
	writeMat(con = paste0(outputPath, "/", file, ".mat"), X = expression, ptime = pseudotime) 
	writeMat(con = paste0(outputPath, "/", file, "_names.mat"), gene_list = rownames(expression))
}
