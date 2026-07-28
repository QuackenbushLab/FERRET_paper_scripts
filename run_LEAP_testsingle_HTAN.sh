#!/bin/bash

# File paths
EXPRESSION_DIR=""
PSEUDOTIME_DIR=""
OUTPUT_DIR=""
EXPRESSION_LARGEST=$EXPRESSION_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
PSEUDOTIME_LARGEST=$PSEUDOTIME_DIR"/pseudotime_RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"
OUTPUT=$OUTPUT_DIR+"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"
echo $(date +%F_%T)
Rscript run_LEAP_single.R  $EXPRESSION_LARGEST $PSEUDOTIME_LARGEST $OUTPUT
echo $(date +%F_%T)
