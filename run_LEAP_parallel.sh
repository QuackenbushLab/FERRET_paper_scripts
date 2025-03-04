#!/bin/bash

# File paths
echo hello
EXPRESSION_DIR="../CPTAC_single_cell_splits_subset_pca_t/"
PSEUDOTIME_DIR="../CPTAC_single_cell_splits_subset_pseudotime/"
OUTPUT_DIR="../LEAP_CPTAC/"
LOGGING_DIR="../LEAP_CPTAC_logs/"
echo $LOGGING_DIR
# For each file in the expression directory, run Panda and save results.
# Create a function that runs the next process only when memory is available.
check_memory_and_run() {
        available_memory=$(free -m | awk '/^Mem:/{print $4}')
        if ((available_memory > 25000)); then
                echo "running LEAP"
		Rscript run_LEAP_single.R  $EXPRESSION_DIR/$1 $PSEUDOTIME_DIR/$2 $OUTPUT_DIR/$1 &> $LOGGING_DIR/$1.txt &
                return 0
        else
                echo "Insufficient memory available. Waiting..."
                return 1
        fi
}

# For each file, check if memory is available, then run.
EXPRESSION_FILES=`ls $EXPRESSION_DIR/C3* | xargs -n 1 basename`
for EFILE in $EXPRESSION_FILES
do
        echo $EFILE
	PSEUDOTIME_FILE=pseudotime_$(echo "${EFILE/%_pca.tsv/.csv}")
        echo $PSEUDOTIME_FILE
        while ! check_memory_and_run $EFILE $PSEUDOTIME_FILE; do
       	        sleep 60
        done
        sleep 60
done

# Wait until all jobs have completed.
wait
echo "All networks ready!"
