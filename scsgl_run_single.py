# Import all necessary packages.
from pysrc.graphlearning import learn_signed_graph
from pysrc.evaluation import auc
import os
import pandas as pd
import numpy as np

# Read expression data.
expression_df = pd.read_csv(sys.argv[1], sep = " ", index_col = 0).transpose()

# Learn graph.
G = learn_signed_graph(expression_df.to_numpy(), pos_density=0.45, neg_density=0.45, 
                                assoc="correlation", gene_names=np.array(expression_df.index))
G.to_csv(sys.argv[2] + ".csv")