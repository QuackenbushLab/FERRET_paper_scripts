#!/bin/bash

# File paths
EXPRESSION_DIR=""
OUTPUT_DIR=""
EXPRESSION_LARGEST=$EXPRESSION_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
OUTPUT=$OUTPUT_DIR"/HTAN_SCENIC/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1"
echo $(date +%F_%T)
sudo Rscript run_scenic.R $EXPRESSION_LARGEST $OUTPUT
echo $(date +%F_%T)

