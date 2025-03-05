#!/bin/bash

# File paths
EXPRESSION_DIR=""
OUTPUT_DIR=""
EXPRESSION_LARGEST=$EXPRESSION_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
OUTPUT=$OUTPUT_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"
mkdir $OUTPUT_DIR

echo $(date +%F_%T)
python3 scgenerai_run_single_tsv.py $EXPRESSION_LARGEST $OUTPUT
echo $(date +%F_%T)
