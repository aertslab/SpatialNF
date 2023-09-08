#!/usr/bin/env python3
import argparse
import scanpy as sc
import squidpy as sq
import os
import pickle

parser = argparse.ArgumentParser(description='''Ligand-receptor interaction analysis.''')

parser.add_argument(
    "h5ad",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "-k", "--cluster_key",
    type=str,
    action="store",
    dest="cluster_key",
    default="louvain",
    help="Key in anndata.AnnData.obs where clustering is stored."
)

parser.add_argument(
    "-t", "--threshold",
    type=float,
    action="store",
    dest="th",
    default=None,
    help="Do not perform permutation test if any of the interacting components is being expressed in less than threshold percent of cells within a given cluster."
)

parser.add_argument(
    "-p", "--complex_policy",
    type=str,
    action="store",
    dest="complex_policy",
    default='min',
    help="Policy on how to handle complexes. 'min' - select gene with the minimum average expression. 'all' - select all possible combinations between ‘source’ and ‘target’ complexes."
)

parser.add_argument(
    "-cm", "--corr_method",
    type=str,
    action="store",
    dest="corr_method",
    default=None,
    help="Correction method for multiple testing. See statsmodels.stats.multitest.multipletests() for valid options."
)

parser.add_argument(
    "-cx", "--corr_axis",
    type=str,
    action="store",
    dest="corr_axis",
    default='clusters',
    help="Axis over which to perform the FDR correction ['interactions', 'clusters']."
)

parser.add_argument(
    "-ca", "--corr_alpha",
    type=float,
    action="store",
    dest="corr_alpha",
    default=0.05,
    help="Significance level for FDR correction."
)

parser.add_argument(
    "-s", "--seed",
    type=int,
    action="store",
    dest="seed",
    default=None,
    help="Random seed for reproducibility."
)

parser.add_argument(
    "-np", "--n_perms",
    type=int,
    action="store",
    dest="n_perms",
    default=1000,
    help="Number of permutations for the permutation test."
)

parser.add_argument(
    "-j", "--n_jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=None,
    help="Number of parallel jobs."
)

parser.add_argument(
    "--not_use_raw",
    choices=['true', 'false'],
    action="store",
    dest="not_use_raw",
    default='false',
    help="Use .raw entry in anndata."
)


parser.add_argument(
    "output_spatial",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_SPATIAL_H5AD = args.h5ad
FILE_PATH_OUT_SPATIAL_BASENAME = os.path.splitext(args.output_spatial.name)[0]

if args.not_use_raw == 'true':
    use_raw_data = False
else:
    use_raw_data = True


# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_SPATIAL_H5AD.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SPATIAL)[0]))

sq.gr.ligrec( adata,
              args.cluster_key,
              threshold=args.th, 
              corr_method=args.corr_method, 
              corr_axis=args.corr_axis, 
              alpha=args.corr_alpha,
              complex_policy=args.complex_policy,
              n_perms=args.n_perms,
              seed=args.seed,
              n_jobs=args.n_jobs,
              use_raw=use_raw_data
)

# Write output
# Usingin pickle as anndata can't handle multiindex (yet)
# https://github.com/theislab/squidpy/issues/409
pickle.dump(adata, open("{}.h5ad".format(FILE_PATH_OUT_SPATIAL_BASENAME), "wb" ))
#adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_SPATIAL_BASENAME))
