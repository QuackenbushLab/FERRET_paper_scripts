#!/bin/bash

# File paths
EXPRESSION_DIR=""
OUTPUT_DIR=""
LOGGING_DIR=""

# For each file in the expression directory, run Panda and save results.
# Create a function that runs the next process only when memory is available.
check_memory_and_run() {
	available_memory=$(free -m | awk '/^Mem:/{print $4}')
	if ((available_memory > 25000)); then
		echo "running SCORPION"
		sudo Rscript run_scorpion.R $EXPRESSION_DIR/$1 $OUTPUT_DIR/$1 &> $LOGGING_DIR/$1.txt &
		return 0
	else
		echo "Insufficient memory available. Waiting..."
		return 1
	fi
}

# For each file, check if memory is available, then run.
EXPRESSION_FILES=`ls $EXPRESSION_DIR/* | xargs -n 1 basename`
for EFILE in $EXPRESSION_FILES; do
	if [[ $EFILE == *original.tsv ]]; then
		echo $EFILE
        	while ! check_memory_and_run $EFILE; do
                	sleep 60
        	done
        	sleep 300
	fi
done
#	for (( i=0; i<${#EXPRESSION_FILES[@]}; i++ )); do
#	if ! test -f $OUTPUT_DIR/${EXPRESSION_FILES[i]}.csv; then
#		echo ${EXPRESSION_FILES[i]}
#		while ! check_memory_and_run ${EXPRESSION_FILES[i]}; do
#			sleep 60
#		done
#sleep 60
#	fi
#done

# Wait until all jobs have completed.
wait
echo "All networks ready!"
