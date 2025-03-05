#!/bin/bash

# File paths
EXPRESSION_DIR=""
PSEUDOTIME_DIR=""
OUTPUT_DIR=""
GRISLI_DIR=""
EXPRESSION_LARGEST=$EXPRESSION_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.tsv"
PSEUDOTIME_LARGEST=$PSEUDOTIME_DIR"/pseudotime_RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"
OUTPUT=$OUTPUT_DIR"/RU682_NORMAL_dense.csv_Macrophages_logcounts.csv_20_1.csv"
cd $GRISLI_DIR
echo $(date +%F_%T)
/usr/share/matlab/bin/matlab -nodisplay -nosplash -nodesktop -r "run_grisli_single_htan('$EXPRESSION_LARGEST', '$PSEUDOTIME_LARGEST', '$OUTPUT'); exit;"
echo $(date +%F_%T)
