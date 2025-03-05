#!/bin/bash

# File paths
EXPRESSION_DIR=""
OUTPUT_DIR=""
EXPRESSION_LARGEST=$EXPRESSION_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
OUTPUT=$OUTPUT_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"

echo $(date +%F_%T)
for CUTOFF in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	julia PIDC_run_single.jl  $EXPRESSION_LARGEST $CUTOFF $OUTPUT
	echo $CUTOFF
done
echo $(date +%F_%T)
