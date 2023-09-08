#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings


# options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "-c",
    "--number-patterns",
    type=int,
    default=5,
    dest='num_patterns',
    help='Number of expected spatial patterns.'
)

parser.add_argument(
    "-l",
    "--length-scale",
    type=float,
    default=None,
    dest='l_value',
    help='Expect length scale of spatial patterns.'
)

parser.add_argument(
    '--estimate-l',
    dest='method_l',
    choices=['none', 'mean', 'median'],
    default='median',
    help="Method for estimating l-value use 'None' to use user-defined l-value."
)

parser.add_argument(
    '--adjust-l',
    dest='adj_l',
    type=float,
    default=0.2,
    help="Fraction to add to estimated l"
)

parser.add_argument(
    "--ncpu",
    type=int,
    default=0,
    dest='ncpu',
    help='Number of CPUs used. 0: does not set environment limit for CPUs'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

lval = args.l_value
cval = args.num_patterns
methodl = args.method_l
adjl = args.adj_l
ncpu = args.ncpu

# env limits for CPUs if specied before importing libs
if ncpu > 0:
    for envkey in ['MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'BLIS_NUM_THREADS']:
        os.environ[envkey] = str(ncpu)
        
import scanpy as sc
import pandas as pd
import SpatialDE
import numpy as np
import scipy


# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")


### main
        
# read anndata
adata = sc.read_h5ad(filename=FILE_PATH_IN.name)

# get count matrix for spatialIDE
if isinstance(adata.X, scipy.sparse.csr.csr_matrix):
    counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
else:
    counts = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)

# get spatial coordinates
if 'X_spatial' in adata.obsm:
    coords = pd.DataFrame(adata.obsm['X_spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
elif 'spatial' in adata.obsm:
    coords = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
    adata.obsm['X_spatial'] = np.float32(adata.obsm['spatial'])[:,:2]
elif 'x' in adata.obs and 'y' in adata.obs:
    coords = pd.DataFrame({'x_coord': adata.obs['x'], 'y_coord': adata.obs['y']}, index=adata.obs_names)
    adata.obsm['X_spatial'] = np.float32(coord)
else:
    raise Exception("VSN ERROR: No spatial coordinate entry found in anndata.")

# get significant genes
signif_res = adata.uns['spatialDE'][adata.uns['spatialDE'].is_signif]

# estimate l value if not specified
if lval == None and methodl == 'none':
    raise Exception("VSN ERROR: '-l' should to be specified if 'estimate_l' is set to 'none.")

if methodl != 'none':

    if methodl == 'median':
        lval = np.median(signif_res['l'])
    elif methodl == 'mean':
        lval = np.mean(signif_res['l'])

    lval = lval + adjl * lval

# run aeh
print("Running automatic expression histology ...")
histology_results, patterns = SpatialDE.aeh.spatial_patterns(coords, counts, signif_res, C=cval, l=lval, verbosity=0)
print("Done.")

# convert column names in to str to make it anndata compliant
patterns.columns = patterns.columns.astype(str)

# add results to uns
adata.uns['spatialDE_AEH_histology_results'] = histology_results
adata.uns['spatialDE_AEH_patterns'] = patterns
# add parameters to anndata
adata.uns['spatialDE_parameters']['l'] = lval
adata.uns['spatialDE_parameters']['c'] = cval

# write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
