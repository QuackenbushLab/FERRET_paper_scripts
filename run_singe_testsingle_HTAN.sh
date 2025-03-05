#!/bin/bash

# File paths
EXPRESSION_LARGEST="RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
cd "../SINGE/SINGE"
mkdir $EXPRESSION_LARGEST
echo $(date +%F_%T)
docker run -v $(pwd):/SINGE -w /SINGE agitter/singe:0.5.1 standalone data4/$EXPRESSION_LARGEST.mat  data4/${EXPRESSION_LARGEST}_names.mat $EXPRESSION_LARGEST default_hyperparameters.txt
echo $(date +%F_%T)
