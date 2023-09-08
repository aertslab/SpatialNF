#!/usr/bin/env python3

import argparse
import os
import scanpy as sc

from umap.parametric_umap import ParametricUMAP
import tensorflow as tf

from matplotlib.colors import to_hex

parser = argparse.ArgumentParser(description='Project spage2vec embedding in rgb color space.')

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
    default=10,
    help="Number of nearest neighbours used."
)

parser.add_argument(
    "--n_epochs",
    type=int,
    action="store",
    dest="n_epochs",
    default=500,
    help="Number of epochs."
)

parser.add_argument(
    "--min_dist",
    type=float,
    action="store",
    dest="min_dist",
    default=0.25,
    help="The effective minimum distance between embedded points."
)

parser.add_argument(
    "--gpu_id",
    type=int,
    action="store",
    dest="gpu_id",
    default=0,
    help="GPU id for training."
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
    "-s", "--seed",
    type=int,
    action="store",
    dest="seed",
    default=0,
    help="Use this integer seed for reproducibility."
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

args = parser.parse_args()

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID";
os.environ["CUDA_VISIBLE_DEVICES"]=str(args.gpu_id)

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

keras_fit_kwargs = {"callbacks": [
    tf.keras.callbacks.EarlyStopping(
        monitor='loss',
        min_delta=10**-2,
        patience=1,
        verbose=1,
    )
]}

reducer = ParametricUMAP(
    n_neighbors=args.n_neigh,
    n_components=3,
    n_epochs=args.n_epochs,
    init='spectral',
    min_dist=args.min_dist,
    spread=1,
    random_state=args.seed,
    verbose=True,
    n_jobs=args.n_cpus,
    keras_fit_kwargs = keras_fit_kwargs,
    unique=True,
    n_training_epochs=1,
    loss_report_frequency=1,
)

embedding = reducer.fit_transform(adata.obsm["spage2vec"])

Y_umap = embedding
Y_umap -= np.min(Y_umap, axis=0)
Y_umap /= np.max(Y_umap, axis=0)

adata.obsm["spage2vec_umap3d"] = Y_umap

adata.obs["spage2vec_umap3d_hex"] = [to_hex(c) for c in adata.obsm['spage2vec_umap3d'].tolist()]

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

