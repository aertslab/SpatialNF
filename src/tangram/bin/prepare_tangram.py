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

### functions

# function for selection marker genes per cluster
def get_best_markerGenes_for_cluster(ad, cluster, max_genes=100, thr_adjp=0.01):

    mgenes = list()
    ngenes = 0
    
    for pos in range(len(ad.uns['rank_genes_groups']['pvals_adj'][cluster])):
        p = ad.uns['rank_genes_groups']['pvals_adj'][cluster][pos]
        if p < thr_adjp:
            mgenes.append(str(ad.uns['rank_genes_groups']['names'][cluster][pos]))
            ngenes += 1
            
            if ngenes==max_genes:
                break
         
    return mgenes


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
    help='Output prefix for spatial and reference h5ad file.'
)


parser.add_argument(
    '-r',
    '--reference',
    dest='ref',
    type=argparse.FileType('r'),
    help='Single cell reference h5ad file containing cell type annotation'
)


parser.add_argument(
    '-n',
    '--number-genes-per-celltype',
    dest='ngenes',
    type=int,
    default=100,
    help="Number of genes per cell type to pick if selected from marker genes (default: '%(default)s')"
)

parser.add_argument(
    '-q',
    '--min-qvalue-celltype',
    dest='qvalue',
    type=float,
    default=0.01,
    help="Min q-value to pick if selected from marker genes (default: %(default)s)"
)

parser.add_argument(
    '-m',
    '--method-gene-selection',
    dest='method',
    choices=['marker_genes', 'list', 'all'],
    default='marker_genes',
    help="Method for selecting genes for cell type mapping (default: '%(default)s')"
)


parser.add_argument(
    '-a',
    '--annotation-celltype',
    dest='anno',
    type=str,
    default='cell.type',
    help="Annotation for selecting from marker genes (default: '%(default)s')"
)


    dest='mode',
    choices=['cells', 'clusters', 'constrained'],
    default='cells',
    help="Mode for cell type mapping (default: '%(default)s')"
)

parser.add_argument(
    "--exp-ref",
    dest='do_exp_ref',
    action="store_true",
    help="Exponentiate reference expression data in case it is log (optional)"
)

parser.add_argument(
    "--exp-spatial",
    dest='do_exp_spatial',
    action="store_true",
    help="Exponentiate spatial expression data in case it is log (optional)"
)

parser.add_argument(
    '--rank-gene-method',
    dest='method_rank_genes',
    choices=['logreg', 't-test', 'wilcoxon', 't-test_overestim_var'],
    default='wilcoxon',
    help="Method for computing 'rank_genes_groups' in reference data if not present (default: '%(default)s')"
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

if args.outmapf:
    FILE_PATH_MAPPING_OUT_BASENAME = os.path.splitext(args.outmapf.name)[0]

### main

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# read scRNAseq file
# Expects h5ad file
try:
    adata_ref = sc.read_h5ad(filename=args.ref.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files for refernce single cell.")

# remove cells and genes with 0 counts everywhere in reference
sc.pp.filter_cells(adata_ref, min_genes=1)
sc.pp.filter_genes(adata_ref, min_cells=1)

# expontiate reference data if log
if args.do_exp_ref:
    adata_ref.X = np.expm1(adata_ref.X)

    
# get marker genes list
gene_list = []

if args.method == 'marker_genes':

    if args.anno not in adata_ref.obs.keys():
        raise Exception("VSN ERROR: Annotation '{}' not found in reference data set.".format(args.anno))
    else:
        if not (hasattr(adata_ref, 'uns') and 'rank_genes_groups' in adata_ref.uns.keys()):
            # compute rank genes groups
            
            # get log for ranking genes
            sc.pp.log1p(adata_ref)

            print("Computing 'rank_genes_groups' ...")
            sc.tl.rank_genes_groups(adata_ref, args.anno, method=args.method_rank_genes)

            # exponentiate for further processing
            adata_ref.X = np.expm1(adata_ref.X)
           
    for cluster in np.unique(adata_ref.obs[args.anno]):
        list2 = get_best_markerGenes_for_cluster(adata_ref, cluster, max_genes=args.ngenes, thr_adjp=args.qvalue)
        gene_list  = gene_list + list(set(list2) - set(gene_list))
                
    if len(gene_list) < 1:
        raise Exception("VSN ERROR: No marker genes selected, adjust selection criteria.")
        
elif args.method == 'all':
    gene_list = adata.var_names

# get intersction of genes
gene_list = list(set(gene_list) & set(adata.var_names) & set(adata_ref.var_names))

if len(gene_list) < 1:
    raise Exception("VSN ERROR: No marker genes selected, adjust selection criteria.")


# tangram specific modification
# add coordinates as obs.x and obs.y
adata_spatial = adata.copy()
adata_spatial.obs['x'] = np.asarray(adata.obsm['spatial'][:,0])
adata_spatial.obs['y'] = np.asarray(adata.obsm['spatial'][:,1])

# exponentiate expression data
if args.do_exp_spatial:
    adata_spatial.X = np.expm1(adata_spatial.X)


# prepare for tangram mapping
tg.pp_adatas(adata_ref, adata_spatial, genes=gene_list)

# run tangram
comput_device = str(args.device)
if comput_device != "cpu":
    comput_device = "cuda: " + comput_device

adata_map = tg.map_cells_to_space(adata_ref, adata_spatial,
                                  device=comput_device,
                                  mode=args.mode,
                                  cluster_label=args.anno,
                                  )

# get mappings
tg.ut.project_cell_annotations(adata_map, adata_spatial, annotation=args.anno)
df_annotations = adata_spatial.obsm["tangram_ct_pred"]

# asssign best cell type
idx = [np.array(df_annotations[df_annotations.index == cellid]).argmax() for cellid in df_annotations.index ]
best_celltypes = df_annotations.columns[idx]
best_celltype_column_name = "tangram_best_" + args.anno


# add cell type annotations use 'n_tangram_' prefix
prefix = "n_tangram_"
df_annotations = df_annotations.add_prefix(prefix)
assert adata.obs.index.equals(df_annotations.index)

# add to obs
df_annotations.index.name = None
adata.obs = pd.concat([adata.obs, df_annotations], axis=1)
adata.obs[best_celltype_column_name] = best_celltypes

# write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

if args.outmapf:
    adata_map.write_h5ad("{}.h5ad".format(FILE_PATH_MAPPING_OUT_BASENAME))
