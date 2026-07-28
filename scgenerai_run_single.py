# Import all necessary packages.
from scGeneRAI import scGeneRAI
import pandas as pd
import sys

# Read expression data.
expression_df = pd.read_csv(sys.argv[1], sep = "\t", index_col = 0).transpose()

# Initialize model and fit data.
model = scGeneRAI()
model.fit(data = expression_df, nepochs = 100, model_depth = 2)

# Predict networks and save.
model.predict_networks(data = expression_df, PATH = sys.argv[2]) 