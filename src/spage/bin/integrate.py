#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
print(sc.__version__)
from SpaGE.principal_vectors import PVComputation
import scipy.stats as st
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='''sn/scRNA-seq and spatial data integration with SpaGE (https://github.com/tabdelaal/SpaGE)''')

parser.add_argument(
    "spatial",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "sc",
    type=argparse.FileType('r'),
    help='Input sn/scRNA-seq reference file.'
)

parser.add_argument(
    "-p", "--npvs",
    type=int,
    action="store",
    dest="n_pv",
    default=30,
    help="Number of principal vectors."
)

parser.add_argument(
    "-t", "--th",
    type=float,
    action="store",
    dest="th",
    default=0.3,
    help="Cosine similarity threshold to remove noisy components."
)

parser.add_argument(
    "-n", "--norm",
    type=bool,
    action="store",
    dest="normalize_sc",
    default=True,
    help="Whether to log-normalize sn/scRNA-seq data before integration."
)


parser.add_argument(
    "output_spatial",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)

parser.add_argument(
    "output_sc",
    type=argparse.FileType('w'),
    help='Output sn/scRNA-seq h5ad file.'
)


args = parser.parse_args()


# Define the arguments properly
FILE_PATH_SPATIAL = args.spatial
FILE_PATH_SC = args.sc
FILE_PATH_OUT_SPATIAL_BASENAME = os.path.splitext(args.output_spatial.name)[0]
FILE_PATH_OUT_SC_BASENAME = os.path.splitext(args.output_sc.name)[0]

# I/O
# Expects h5ad file
try:
    ad_spatial = sc.read_h5ad(filename=FILE_PATH_SPATIAL.name)
    if ad_spatial.X.min() < 0:
        # Data are probably scaled
        ad_spatial.X = ad_spatial.raw.to_adata().X
    print(ad_spatial)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SPATIAL)[0]))

try:
    ad_sc = sc.read_h5ad(filename=FILE_PATH_SC.name)
    print(ad_sc)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SC)[0]))

# Log-normalize datasets
if args.normalize_sc:
    sc.pp.normalize_total(ad_sc)
    sc.pp.log1p(ad_sc) 

#
# Project data in common latent space by domain adaptation
#

n_pv = args.n_pv

# Extract and scale expression matrices
df_spatial = ad_spatial.to_df()
df_sc = ad_sc.to_df()
df_spatial = pd.DataFrame(st.zscore(df_spatial.values), index=df_spatial.index, columns=df_spatial.columns)
df_sc = pd.DataFrame(st.zscore(df_sc.values), index=df_sc.index, columns=df_sc.columns)

# Subset sn/scRNA-seq expression matrix on genes shared with spatial data
Common_data = df_sc[np.intersect1d(df_spatial.columns, df_sc.columns)]
print(Common_data)
# Compute shared latent space
pv_FISH_RNA = PVComputation(
        n_factors = n_pv,
        n_pv = n_pv,
        dim_reduction = 'pca',
        dim_reduction_target = 'pca'
)
pv_FISH_RNA.fit(Common_data.values, df_spatial[Common_data.columns].values)

S = pv_FISH_RNA.source_components_.T

# Remove noisy components based on threshold
Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > args.th)
S = S[:,0:Effective_n_pv]
# Project data into shared latent space
Common_data_t = Common_data.dot(S)
FISH_exp_t = df_spatial[Common_data.columns].dot(S)

ad_spatial.obsm['X_spage'] = FISH_exp_t.values
ad_sc.obsm['X_spage'] = Common_data_t.values
# Write outputs
ad_spatial.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_SPATIAL_BASENAME))
ad_sc.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_SC_BASENAME))
