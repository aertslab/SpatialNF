#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import numpy as np
import pickle

parser = argparse.ArgumentParser(description='''knn-imputation''')

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
    "-k", "--knn",
    type=int,
    action="store",
    dest="k_nn",
    default=15,
    help="Number of k-nearest neighbors for imputation."
)

parser.add_argument(
    "-j", "--n-jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=1,
    help="The number of jobs. When set to -1, automatically uses the number of cores."
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output file.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_SPATIAL = args.spatial
FILE_PATH_SC = args.sc
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

#
# Run k-nn imputation
#

nbrs = NearestNeighbors(n_neighbors=args.k_nn, algorithm='auto',
                            metric = 'cosine', n_jobs=args.n_jobs).fit(ad_sc.obsm['X_spage'])
res = nbrs.kneighbors(ad_spatial.obsm['X_spage'])

# Write outputs
pickle.dump( res, open( "{}.pkl".format(FILE_PATH_OUT_BASENAME), "wb" ) )
