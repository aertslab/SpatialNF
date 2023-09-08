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
import random

### functions

# function for selection marker genes per cluster
def get_best_markerGenes_for_cluster(ad, cluster, sparsity_key, max_genes=100, thr_adjp=0.01, max_sparse=0.5):

    mgenes = list()
    ngenes = 0
    
    for pos in range(len(ad.uns['rank_genes_groups']['pvals_adj'][cluster])):

        p = ad.uns['rank_genes_groups']['pvals_adj'][cluster][pos]
        gene_name = str(ad.uns['rank_genes_groups']['names'][cluster][pos])
        idx_gene = ad.uns['var_getIndex'][gene_name]
        
        sparsity = adata_ref.uns[sparsity_key][cluster][idx_gene]
        
        if p < thr_adjp and sparsity <= max_sparse:
            mgenes.append(gene_name)
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
    help='Output spatial h5ad file.'
)


parser.add_argument(
    "--output-mapping",
    type=argparse.FileType('w'),
    dest='outmapf',
    help='Output mapping h5ad file. (Optional)'
)


parser.add_argument(
    '-r',
    '--reference',
    dest='ref',
    type=argparse.FileType('r'),
    help='Single cell reference h5ad file containing cell type annotation'
)

parser.add_argument(
    '-l',
    '--list_genes',
    dest='list_genes',
    default=None,
    type=argparse.FileType('r'),
    help='File containing marker genes, one gene per row. To use gene list, set --method-gene-selection=list (optional)'
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
    help="Min q-value to pick gene if selected from marker genes (default: %(default)s)"
)

parser.add_argument(
    '-s',
    '--max-sparsity',
    dest='maxsparse',
    type=float,
    default=0.5,
    help="Max sparsity per cluster to pick gene if selected from marker genes (default: %(default)s)"
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
    '-d',
    '--device',
    default='0',
    help="Device used for computing cell type mapping, e.g. 'cpu' for CPU, '0' for 'cuda: 0', ... for GPU (default: '%(default)s')" 
)

parser.add_argument(
    '-a',
    '--annotation-celltype',
    dest='anno',
    type=str,
    default='cell.type',
    help="Annotation for selecting from marker genes (default: '%(default)s')"
)


parser.add_argument(
    '--mode',
    dest='mode',
    choices=['cells', 'clusters', 'constrained'],
    default='cells',
    help="Mode for cell type mapping (default: '%(default)s')"
)

parser.add_argument(
    '--seed',
    dest='seed',
    type=int,
    default=None,
    help="Initialize RNG to seed (default: '%(default)s')"
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

# initialize RNG if specified
if args.seed:
    random.seed(args.seed)

# get marker genes list
gene_list = []

if args.method == 'marker_genes':

    if args.anno not in adata_ref.obs.keys():
        raise Exception("VSN ERROR: Annotation '{}' not found in reference data set.".format(args.anno))
    else:
        if not (hasattr(adata_ref, 'uns') and 'rank_genes_groups' in adata_ref.uns.keys()):
            raise Exception("VSN ERROR: 'rank_genes_groups' not found in reference data set. Should have been computed by workflow.")
    if type(adata_ref.uns['rank_genes_groups']['pvals_adj']) == dict:
        cluster_names = adata_ref.uns['rank_genes_groups']['pvals_adj'].keys()
    elif type(adata_ref.uns['rank_genes_groups']['pvals_adj']) == np.ndarray:
        cluster_names = adata_ref.uns['rank_genes_groups']['pvals_adj'].dtype.names
    else:
        raise Exception("VSN ERROR: .uns['rank_genes_groups']['pvals_adj'] in reference has to be of type 'dict' or 'numpy.ndarray' including 'dtype.names'")

    if 'var_getIndex' not in adata_ref.uns.keys():
        raise Exception("VSN ERROR: .uns['var_getIndex'] not in reference")
    sparse_key = args.anno + '_sparsity'
    if sparse_key not in adata_ref.uns.keys():
        raise Exception("VSN ERROR: .uns['" + sparse_key + "'] not in reference")
    
    for cluster in cluster_names:
        list2 = get_best_markerGenes_for_cluster(adata_ref, cluster, sparsity_key=sparse_key, max_genes=args.ngenes, thr_adjp=args.qvalue, max_sparse=args.maxsparse)
        gene_list  = gene_list + list(set(list2) - set(gene_list))
                
    if len(gene_list) < 1:
        raise Exception("VSN ERROR: No marker genes selected, adjust selection criteria.")
        
elif args.method == 'all':
    gene_list = adata_spatial.var_names
    
elif args.method == 'list':
    if args.list_genes is not None:
        try:
            with open(args.list_genes.name, 'r') as f:
                lines = f.readlines()
                gene_list = [ line.strip() for line in lines ]
        except IOError:
            raise Exception("VSN ERROR: Failed reading marker gene file {}.".format(args.list_genes))
    else:
        raise Exception("VSN ERROR: No marker gene file specified. '--list_genes' has to be specified if --method-gene-selection=list.")

# get intersection of genes
gene_list = list(set(gene_list) & set(adata_spatial.var_names) & set(adata_ref.var_names))

if len(gene_list) < 1:
    raise Exception("VSN ERROR: No marker genes selected, adjust selection criteria.")

# prepare for tangram mapping
tg.pp_adatas(adata_ref, adata_spatial, genes=gene_list)

# run tangram
comput_device = str(args.device)
if comput_device == "any":
    comput_device = "cuda: " + str(torch.cuda.current_device())
elif comput_device != "cpu":
    comput_device = "cuda: " + comput_device

adata_map = tg.map_cells_to_space(adata_ref, adata_spatial,
                                  device=comput_device,
                                  mode=args.mode,
                                  cluster_label=args.anno,
                                  )

# write output
adata_map.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
