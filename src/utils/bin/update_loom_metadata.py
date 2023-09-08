#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import scanpy as sc
import numpy as np
from loomxpy._loomx import LoomX
from loomxpy._io._read import read_loom

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input_loom",
    type=argparse.FileType('r'),
    help='The path to the input loom file.'
)

parser.add_argument(
    "input_h5ad",
    type=argparse.FileType('r'),
    help='The path to the input h5ad file storing the additional metadata to add to the given input loom file.'
)


parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output loom file.'
)

parser.add_argument(
    '-a', '--axis',
    type=str,
    dest="axis",
    required=True,
    default='observation',
    help='The axis to which the metadatas will be added to. Possible values: [feature,observation].'
)

parser.add_argument(
    '-m', "--metric-keys",
    type=str,
    dest="metric_keys",
    required=False,
    help='Comma-separated list of metric metadata to add.'
)

parser.add_argument(
    '-t', '--annotation-keys',
    type=str,
    dest="annotation_keys",
    required=False,
    help='Comma-separated list of annotation metadata to add.'
)

parser.add_argument(
    '-e', '--embedding-keys',
    type=str,
    dest="embedding_keys",
    required=False,
    help="Comma-separated list of embedding metadata to add."
)

parser.add_argument(
    '-c', '--clustering-keys',
    type=str,
    dest="clustering_keys",
    help="Comma-separated list of clustering metadata to add."
)


args = parser.parse_args()

FILE_PATH_IN_LOOM = args.input_loom.name
FILE_PATH_IN_H5AD = args.input_h5ad.name
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    loom = read_loom(file_path=FILE_PATH_IN_LOOM,
        force_conversion={"annotations": True}
    )
except IOError:
    raise Exception("VSN ERROR: Can only handle scope loom files.")
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN_H5AD)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

#
# Update loom metadata.
#
if args.axis == 'observation':
    # Add metric metadata
    if args.metric_keys:
        metric_keys = args.metric_keys.split(',')
        for k in metric_keys:
            loom.modes.rna.o.metrics[k] = pd.DataFrame({k: adata.obs.loc[:,k]}, index=adata.obs.index)

    # Add annotation metadata
    if args.annotation_keys:
        annotation_keys = args.annotation_keys.split(',')
        for k in annotation_keys:
            loom.modes.rna.o.annotations[k] = pd.DataFrame({k: adata.obs.loc[:,k]}, index=adata.obs.index)

    # Add embedding metadata
    if args.embedding_keys:
        embedding_keys = args.embedding_keys.split(',')
        for k in embedding_keys:
            loom.modes.rna.o.embeddings[k] = pd.DataFrame({k: adata.obsm[k]}, index=adata.obs.index)

    # Add clustering metadata
    if args.clustering_keys:
        # Add clustering to loom not working propely in LoomXpy
        raise Exception("VSN ERROR: Updating the clustering metadata is currently not implemented.")
#        clustering_keys = args.clustering_keys.split(',')
#        for k in clustering_keys:
#            if adata.obs.loc[:,k].dtype == int or float:
#                loom.modes.rna.o.clusterings[k] = pd.DataFrame({k: adata.obs.loc[:,k]}, index=adata.obs.index)
#            else:
#                cluster_labels = np.unique(adata.obs.loc[:,k])
#                cl2label = dict(zip(cluster_labels, np.arange(cluster_labels.shape[0])))
#                loom.modes.rna.o.clusterings[k] = pd.DataFrame({k: adata.obs.loc[:,k].map(cl2label)}, index=adata.obs.index)

elif args.axis == 'feature':
    # Add metric metadata
    if args.metric_keys:
        metric_keys = args.metric_keys.split(',')
        for k in metric_keys:
            metric_df = pd.DataFrame(index=adata.raw.var.index).merge(adata.var[k], left_index=True, right_index=True, how="left")
            loom.modes.rna.f.metrics[k] = metric_df

    # Add annotation metadata
    if args.annotation_keys:
        annotation_keys = args.annotation_keys.split(',')
        for k in annotation_keys:
            annotation_df = pd.DataFrame(index=adata.raw.var.index).merge(adata.var[k], left_index=True, right_index=True, how="left").fillna(value=False)
            loom.modes.rna.f.annotations[k] = annotation_df
else:
    raise Exception(f"Cannot update the {args.axis}-based metadata.")


# I/O
loom.modes.rna.export(filename="{}.loom".format(FILE_PATH_OUT_BASENAME), output_format="scope_v1")
