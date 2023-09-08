#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import pickle

parser = argparse.ArgumentParser(description='''Label tranfer of sc/snRNA-seq annotations to spatial data''')

parser.add_argument(
    "spatial",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "sc",
    type=argparse.FileType('r'),
    help='Input sn/scRNA-seq h5ad file.'
)

parser.add_argument(
    "knn_input",
    type=argparse.FileType('r'),
    help='knn-imputation result file.'
)

parser.add_argument(
    "-k", "--label-key",
    type=str,
    action="store",
    dest="key",
    default='cell_type',
    help="adata.obs column name containing labels to transfer."
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file with label transfered.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_SPATIAL = args.spatial
FILE_PATH_SC = args.sc
FILE_PATH_KNN_INPUT = args.knn_input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    ad_spatial = sc.read_h5ad(filename=FILE_PATH_SPATIAL.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SPATIAL)[0]))

try:
    ad_sc = sc.read_h5ad(filename=FILE_PATH_SC.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SC)[0]))

try:
    knn_res = pickle.load( open( FILE_PATH_KNN_INPUT.name, "rb" ) )
except IOError:
    raise Exception("Wrong input format. Expects binary pickle files, got .{}".format(os.path.splitext(FILE_PATH_KNN_INPUT)[0]))

#
# Run label transfering
#
distances, indices = knn_res
labels = ad_sc.obs.loc[:, args.key].unique()
ct_df = pd.DataFrame(np.zeros((ad_spatial.shape[0], len(labels))), index=ad_spatial.obs.index, columns=labels)
for j in range(0, ad_spatial.shape[0]):
    weights = 1 - (distances[j,:][distances[j,:]<1]) / (np.sum(distances[j,:][distances[j,:]<1]))
    weights = weights / (len(weights)-1)
    for c in ad_sc.obs.iloc[indices[j,:], ad_sc.obs.columns.get_loc(args.key)].unique():
        c_idx = np.where(ad_sc.obs.iloc[indices[j,:], ad_sc.obs.columns.get_loc(args.key)].values==c)[0]
        ct_df.iloc[j, ct_df.columns.get_loc(c)] = np.sum(weights[c_idx])

ad_spatial.obs[args.key] = ct_df.idxmax(axis=1)
ad_spatial.obs[args.key + '_dist'] = ct_df.max(axis=1)

# Write outputs
ad_spatial.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
