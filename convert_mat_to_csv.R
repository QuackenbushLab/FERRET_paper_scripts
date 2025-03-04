if(!require("R.matlab")){
	install.packages("R.matlab")
}
library("R.matlab")

matPath = "/home/ubuntu/SINGE/SINGE"
outputPath = "/home/ubuntu/SINGE_csv"
dir.create(outputPath)

# Loop through each file and process it
for(dir in list.dirs(matPath)){

	# The ID and lambda values change. For a given ID, the replicate number is constant.
	# The file name contains the name of the same without the "_tsv".
	result = readMat(con = paste0(matPath, "/", dir, "/", "AdjMatrix_data1_C3N-01904-02_Macrophages_1_10_pcaptsvpmat_ID_9_lambda_0_replicate_1.mat"))
	
	# Save as CSV.
	str(result)
	write.csv(result, con = paste0(outputPath, "/", dir, ".csv"))
}
