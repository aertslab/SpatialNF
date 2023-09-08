#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import scanpy as sc
import pandas as pd
import numpy as np


# options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "final",
    type=argparse.FileType('r'),
    help='Proccesed h5ad file patterns should be added to.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_FINAL = args.final
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

### main

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# Expects h5ad file
try:
    adata_final = sc.read_h5ad(filename=FILE_PATH_FINAL.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# to prevent pot. encoding issue substitute tuple entries with strings
for obskey in adata_final.obs.keys():
    if isinstance(adata_final.obs[obskey][0], tuple):
        adata_final.obs[obskey] = [ str(tupl) for tupl in adata_final.obs[obskey] ]
for varkey in adata_final.var.keys():
    if isinstance(adata_final.var[varkey][0], tuple):
        adata_final.var[varkey] = [ str(tupl) for tupl in adata_final.var[varkey] ]

# add patterns as obs variables to used as SCope metric
for pattern in adata.uns['spatialDE_AEH_patterns'].keys():
    patname = 'n_spatialDE_AEH_pattern' + pattern
    adata_final.obs[patname] = [ np.exp(v) for v in adata.uns['spatialDE_AEH_patterns'][pattern] ]

# write output
adata_final.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
