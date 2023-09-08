#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import squidpy as sq

parser = argparse.ArgumentParser(description='''Compute neighborhood enrichment analysis by permutation test.''')

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
    "-c", "--connectivity_key",
    type=str,
    action="store",
    dest="connectivity_key",
    default="spatial",
    help="Key in anndata.AnnData.obsp where spatial connectivities are stored."
)

parser.add_argument(
    "-p", "--n_perms",
    type=int,
    action="store",
    dest="n_perms",
    default=1000,
    help="Number of permutations for the permutation test."
)

parser.add_argument(
    "-m", "--ripley_mod",
    type=str,
    action="store",
    dest="ripley_mod",
    default='L',
    help="Which Ripley’s statistic to compute."
)


parser.add_argument(
    "-n", "--n_steps",
    type=int,
    action="store",
    dest="n_steps",
    default=50,
    help="Number of distance thresholds at which cluster co-occurrence is computed."
)

parser.add_argument(
    "-nb", "--n_neigh",
    type=int,
    action="store",
    dest="n_neigh",
    default=2,
    help="Number of neighbors to consider for the KNN graph."
)

parser.add_argument(
    "-q", "--n_simulations",
    type=int,
    action="store",
    dest="n_simulations",
    default=100,
    help="How many simulations to run for computing p-values for Ripley's statistics."
)

parser.add_argument(
    "-o", "--n_observations",
    type=int,
    action="store",
    dest="n_observations",
    default=1000,
    help="How many observations to generate for the Spatial Poisson Point Process in Ripley's statistics."
)

parser.add_argument(
    "-d", "--max_dist",
    type=float,
    action="store",
    dest="max_dist",
    default=None,
    help="Maximum distances for the support in Ripley's statistics."
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
    "-j", "--n_jobs",
    type=int,
    action="store",
    dest="n_jobs",
    default=None,
    help="Number of parallel jobs."
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

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_SPATIAL_H5AD.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SPATIAL)[0]))

# neighborhood enrichment
sq.gr.nhood_enrichment( adata,
    cluster_key = args.cluster_key,
    connectivity_key = args.connectivity_key,
    n_perms = args.n_perms,
    seed = args.seed,
    n_jobs = args.n_jobs)

# interaction matrix
sq.gr.interaction_matrix( adata,
    cluster_key = args.cluster_key,
    connectivity_key = args.connectivity_key)

# cluster co-occurence probabilities
sq.gr.co_occurrence( adata,
    cluster_key = args.cluster_key,
    spatial_key = args.connectivity_key,
    n_steps = args.n_steps,
    n_jobs = args.n_jobs)

# centrality scores
sq.gr.centrality_scores( adata,
    cluster_key = args.cluster_key,
    connectivity_key = args.connectivity_key,
    n_jobs = args.n_jobs)

# Ripley’s statistics
sq.gr.ripley( adata,
    cluster_key = args.cluster_key,
    mode = args.ripley_mod,
    spatial_key = args.connectivity_key,
    n_steps = args.n_steps,
    n_neigh = args.n_neigh,
    n_simulations = args.n_simulations,
    n_observations = args.n_observations,
    max_dist = args.max_dist,
    seed = args.seed,
)

# Write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_SPATIAL_BASENAME))
