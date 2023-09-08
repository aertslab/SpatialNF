#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import squidpy as sq
import anndata as ad

parser = argparse.ArgumentParser(description='''Compute spatial neighbors graph.''')

parser.add_argument(
    "h5ad",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "-k", "--spatial_key",
    type=str,
    action="store",
    dest="spatial_key",
    default="spatial",
    help="Key in anndata.AnnData.obsm where spatial coordinates are stored."
)

parser.add_argument(
    "-t", "--coord_type",
    type=str,
    action="store",
    dest="coord_type",
    default=None,
    help="Type of coordinate system (valid options: ['grid','generic',None]. None - ‘grid’ if spatial_key is in anndata.AnnData.uns with n_neighs = 6 (Visium), otherwise use ‘generic’."
)

parser.add_argument(
    "-n", "--n_neighs",
    type=int,
    action="store",
    dest="n_neighs",
    default=6,
    help="Number of neighborhoods."
)

parser.add_argument(
    "-r", "--radius",
    type=float,
    action="store",
    dest="radius",
    default=None,
    help="Compute the graph based on neighborhood radius."
)

parser.add_argument(
    "-d", "--delaunay",
    type=bool,
    action="store",
    dest="delaunay",
    default=False,
    help="Whether to compute the graph from Delaunay triangulation. Only used when coord_type = 'generic'."
)

parser.add_argument(
    "-s", "--n_rings",
    type=int,
    action="store",
    dest="n_rings",
    default=1,
    help="Number of rings of neighbors for grid data. Only used when coord_type = 'grid'."
)

parser.add_argument(
    "-l", "--set_diag",
    type=bool,
    action="store",
    dest="set_diag",
    default=False,
    help="Whether to set the diagonal of the spatial connectivities to 1.0."
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

# If multiple sample generate disjoint spatial graphs

sample_ids = adata.obs.sample_id.unique()
if len(sample_ids)>1:
    ad_l = []
    for s in adata.obs.sample_id.unique():
        adata_tmp = adata[adata.obs.sample_id==s].copy()
        sq.gr.spatial_neighbors( adata_tmp,
                         spatial_key=args.spatial_key, 
                         coord_type=args.coord_type,
                         n_neighs=args.n_neighs, 
                         radius=args.radius,
                         delaunay=args.delaunay,
                         n_rings=args.n_rings,
                         set_diag=args.set_diag,
                         key_added=args.spatial_key,
                         copy=False )
        ad_l.append(adata_tmp)
    adata = ad.concat(ad_l, uns_merge='same', pairwise=True)
else:
    sq.gr.spatial_neighbors( adata,
                         spatial_key=args.spatial_key, 
                         coord_type=args.coord_type,
                         n_neighs=args.n_neighs, 
                         radius=args.radius,
                         delaunay=args.delaunay,
                         n_rings=args.n_rings,
                         set_diag=args.set_diag,
                         key_added=args.spatial_key,
                         copy=False )


# Write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_SPATIAL_BASENAME))
