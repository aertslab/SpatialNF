#!/usr/bin/env python3

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.spatial import cKDTree as KDTree

parser = argparse.ArgumentParser(description='Boostrap psuedo cellular profiles from spage2vec embedding.')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file with spatial graph.'
)

parser.add_argument(
    "--n_neigh",
    type=int,
    action="store",
    dest="n_neigh",
    default=100,
    help="Number of neighboring node to bootstrap from."
)

parser.add_argument(
    "--n_cpus",
    type=int,
    action="store",
    dest="n_cpus",
    default=1,
    help="Number of jobs for multiprocessing."
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

args = parser.parse_args()
# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

kdT = KDTree(adata.obsm["spage2vec"])
res = kdT.query(adata.obsm["spage2vec"], k=args.n_neigh, n_jobs=args.n_cpus)

spot_geneID = adata.obs.gene.values
gene_list = adata.obs.gene.unique()

df_pc_exp = pd.DataFrame(np.zeros((len(res[1]), len(gene_list))), columns=gene_list)
col_idx_d = dict(zip(gene_list, np.arange(gene_list.shape[0])))
for spot in range(len(res[1])):
    geneID, counts = np.unique(spot_geneID[res[1][spot]], return_counts=True)
    geneID_idx = [col_idx_d[x] for x in geneID]
    df_pc_exp.iloc[spot, geneID_idx] = counts

adata.varm['pc_var_names'] = np.array(df_pc_exp.columns.values)
adata.obsm['pc'] = df_pc_exp.values

adata_pc = AnnData(adata.obsm["pc"], obs= adata.obs)
adata_pc.var_names = adata.varm["pc_var_names"]
adata_pc.obsm = adata.obsm

# If spatial dataset arrange 'X_spatial' embedding of each sample in a grid for SCope visualization
if hasattr(adata_pc, 'obsm') and 'spatial' in adata.obsm.keys():
    samples = adata_pc.obs.sample_id.unique()
    n_samples = len(samples)
    n_x = np.ceil(np.sqrt(n_samples)) if np.sqrt(n_samples) >= np.floor(np.sqrt(n_samples)) else np.floor(np.sqrt(n_samples))
    n_y = np.ceil(np.sqrt(n_samples)) if np.sqrt(n_samples) >= np.floor(np.sqrt(n_samples))+0.5 else np.floor(np.sqrt(n_samples))

    # get spacing between samples from first sample to define offset on fixed fraction of width/height, e.g. 10%
    frac_offset = 0.1
    s = samples[0]
    s1_width = adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,0].max()-adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,0].min()
    s1_height = adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,1].max()-adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,1].min()
    offset = np.max([s1_width, s1_height]) * frac_offset

    x = 0
    x_max = y_max = 0
    X = []
    Y = []
    for s in samples:
        # Reset origin in (0,0)
        x_s = adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,0]-adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,0].min()
        y_s = adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,1]-adata_pc[adata_pc.obs.sample_id==s].obsm['spatial'][:,1].min()

        X.append(x_s + x_max)
        Y.append(y_s + y_max)

        x_max = x_max + x_s.max() + offset
        x = x + 1
        if x>n_x-1:
            x = 0
            x_max = 0
            y_max = y_max - y_s.max() - offset
    adata_pc.obsm['X_spatial'] = np.array([np.concatenate(X), np.concatenate(Y)]).T

    # center data for SCope
    avg_coords = adata_pc.obsm['X_spatial'].sum(axis=0)/adata_pc.obsm['X_spatial'].shape[0]
    adata_pc.obsm['X_spatial'] = adata_pc.obsm['X_spatial'] - avg_coords
    # scale data for SCope such that X axis is between [-10,10]
    scfactor = 20/(np.max(adata_pc.obsm['X_spatial'][:,0])-np.min(adata_pc.obsm['X_spatial'][:,0]))
    adata_pc.obsm['X_spatial'] = adata_pc.obsm['X_spatial'] * scfactor


#adata_pc.obsm["spage2vec"] = adata.obsm["spage2vec"]
#adata_pc.obsm["spatial"] = adata.obsm["spatial"]
#adata_pc.obs["sample_id"] = adata.obs["sample_id"]

# Log-normalize pseudo-cell expression profiles
#sc.pp.normalize_total(adata_pc)
#sc.pp.log1p(adata_pc)

adata_pc.write("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
