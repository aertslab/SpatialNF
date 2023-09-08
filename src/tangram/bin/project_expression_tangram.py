#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import torch
import tangram as tg
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ann

### options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)


parser.add_argument(
    '--mapping',
    type=argparse.FileType('r'),
    dest='map',
    help='Input mapping h5ad file.'
)


parser.add_argument(
    '-r',
    '--reference',
    dest='ref',
    type=argparse.FileType('r'),
    help='Single cell reference h5ad file containing cell type annotation'
)

parser.add_argument(
    '-a',
    '--annotation-celltype',
    dest='anno',
    type=str,
    default='cell.type',
    help="Annotation for selecting from marker genes (default: '%(default)s')"
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


### main

# I/O
# Expects h5ad file
try:
    adata_spatial = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# read scRNAseq file
# Expects h5ad file
try:
    adata_ref = sc.read_h5ad(filename=args.ref.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files for refernce single cell.")

# read mapping file
# Expects h5ad file
try:
    adata_map = sc.read_h5ad(filename=args.map.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files for tangram mapping.")

# save original gene names as tangram sets them to lowercase
dict_genename = {}
for g in adata_ref.var.index:
    dict_genename[g.lower()] = g

# get mappings
ad_ge = tg.project_genes(adata_map, adata_ref)

# create new anndata
adata_gex = ann.AnnData(X=ad_ge.X, obs=adata_spatial.obs, var=ad_ge.var, uns=adata_spatial.uns, obsm=adata_spatial.obsm, obsp=adata_spatial.obsp)

# get original gene names from scRNAseq as tangram appears to rename them to lowercase etc.
adata_gex.var['tangram_gene'] = adata_gex.var.index

new_index = []
for g in adata_gex.var.index:
    if g in dict_genename:
        new_index.append(dict_genename[g])
    else:
        new_index.append(g)
adata_gex.var.index = new_index
adata_gex.var['Gene'] = new_index

# write output
adata_gex.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
