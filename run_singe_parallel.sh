#!/bin/bash

# File paths
EXPRESSION_DIR="../../CPTAC_single_cell_splits_subset_pca/"
cd "../SINGE/SINGE"
# For each file in the expression directory, run Panda and save results.
# Create a function that runs the next process only when memory is available.
check_memory_and_run() {
	available_memory=$(free -m | awk '/^Mem:/{print $4}')
	if ((available_memory > 25000)); then
		echo "running SINGE"
		docker run -v $(pwd):/SINGE -w /SINGE agitter/singe:0.5.1 standalone data1/$1.mat  data1/$1_names.mat $1 default_hyperparameters.txt &
		return 0
	else
		echo "Insufficient memory available. Waiting..."
		return 1
	fi
}

# For each file, check if memory is available, then run.
EXPRESSION_FILES=`ls $EXPRESSION_DIR/C3* | xargs -n 1 basename`
EFILE_SUBSET_INC=$(comm -23 <(sort <(find "$EXPRESSION_DIR" -type f -printf "%f.txt\n")) <(sort <(find "$LOGGING_DIR" -type f -printf "%f\n")))
echo $EFILE_SUBSET_INC
for EFILE in $EXPRESSION_FILES
#for FILE in $EFILE_SUBSET_INC
do
	#EFILE=$(awk -F ".txt" '{ for(i=1;i<=NF;i++) print $i }' <<< $FILE)
	echo $EFILE
	while ! check_memory_and_run $EFILE; do
		sleep 60
	done
	sleep 300
done

# Wait until all jobs have completed.
wait
echo "All networks ready!"
