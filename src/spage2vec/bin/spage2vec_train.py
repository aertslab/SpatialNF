#!/usr/bin/env python3

import argparse
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import datetime
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData
import anndata as ann

import numpy as np
import networkx as nx
#import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.spatial import cKDTree as KDTree
from scipy.spatial.distance import euclidean

parser = argparse.ArgumentParser(description='Train spage2vec on spatial gene expression graph.')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file with spatial graph.'
)

parser.add_argument(
    "-f", "--feature_key",
    type=str,
    action="store",
    dest="feature_key",
    default="gene",
    help="Key in anndata.AnnData.obs where node features are stored."
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
    "-r", "--radius",
    type=float,
    action="store",
    dest="radius",
    default=None,
    help="Compute the graph based on neighborhood radius. Or specify neigh_perc."
)

parser.add_argument(
    "-p", "--neigh_perc",
    type=int,
    action="store",
    dest="neigh_perc",
    default=97,
    help="Compute the spatial graph such that <neigh_perc> percent of the graph nodes have at least one neighbour."
)

parser.add_argument(
    "-d", "--delaunay",
    type=bool,
    action="store",
    dest="delaunay",
    default=False,
    help="Whether to compute the graph from Delaunay triangulation."
)

parser.add_argument(
    "--min_conn_comp_size",
    type=int,
    action="store",
    dest="min_conn_comp_size",
    default=3,
    help="Size of the smallest connected component allowed."
)

parser.add_argument(
    "--num_of_walks",
    type=int,
    action="store",
    dest="num_of_walks",
    default=1,
    help="Number of random walks for sampling neighbours."
)

parser.add_argument(
    "--walk_length",
    type=int,
    action="store",
    dest="walk_length",
    default=2,
    help="Length of random walks."
)

parser.add_argument(
    "--batch_size",
    type=int,
    action="store",
    dest="batch_size",
    default=32,
    help="Batch size for model training."
)

parser.add_argument(
    "--num_epochs",
    type=int,
    action="store",
    dest="num_epochs",
    default=50,
    help="Number of epochs for model training."
)

parser.add_argument(
    "--k1",
    type=int,
    action="store",
    dest="k1",
    default=10,
    help="Number of nodes to sample in the first hop."
)

parser.add_argument(
    "--k2",
    type=int,
    action="store",
    dest="k2",
    default=5,
    help="Number of nodes to sample in the second hop."
)

parser.add_argument(
    "--l1",
    type=int,
    action="store",
    dest="l1",
    default=20,
    help="Size (num. of neurons) of layer one."
)

parser.add_argument(
    "--l2",
    type=int,
    action="store",
    dest="l2",
    default=20,
    help="Size (num. of neurons) of layer two."
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
    "--n_cpus",
    type=int,
    action="store",
    dest="n_cpus",
    default=1,
    help="Number of jobs for multiprocessing."
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
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "spage2vec_model",
    type=argparse.FileType('w'),
    help='Output spage2vec model file.'
)

args = parser.parse_args()

#os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID";
#os.environ["CUDA_VISIBLE_DEVICES"]=str(args.gpu_id)

import tensorflow as tf
from tensorflow.keras.utils import to_categorical
from tensorflow import keras

from stellargraph import StellarGraph
from stellargraph.data import UnsupervisedSampler
from stellargraph.layer import GraphSAGE, link_classification
from stellargraph.layer.graphsage import AttentionalAggregator
from stellargraph.mapper import GraphSAGELinkGenerator, GraphSAGENodeGenerator

import random as rn


# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
FILE_PATH_MODEL_BASENAME = args.spage2vec_model.name

# I/O
# Expects h5ad file
try:
    df = pd.read_csv(FILE_PATH_IN.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_IN)[0]))

################################################################################
def computeNeighborDistance(ad, th, sample_key='sample_id', spatial_key='spatial'):
    d_list=[]
    for s in ad.obs.loc[:,sample_key].unique():
        kdT = KDTree(ad[ad.obs.loc[:,sample_key]==s].obsm[spatial_key].toarray())
        d,i = kdT.query(ad[ad.obs.loc[:,sample_key]==s].obsm[spatial_key].toarray(),k=2)
        d_list.append(d)
    d = np.vstack(d_list)
    print(d.shape)
    #plt.hist(d[:,1],bins=200);
    #plt.axvline(x=np.percentile(d[:,1],th),c='r')
    print(np.percentile(d[:,1],th))
    d_th = np.percentile(d[:,1],th)
    return d_th

def buildGraph(ad, use_rep='X'):
    ad = ad.copy()
    ad.obs.reset_index(inplace=True, drop=True)
    G = nx.Graph()
    node_removed = []
    
    res = np.array(ad.obsp['spatial_connectivities'].nonzero()).T
    d_res = np.array(ad.obsp['spatial_distances'][res[:,0],res[:,1]])[0]
    res = [(x[0],x[1]) for x in list(res)]

    if use_rep == 'X':
        features = ad.X
    else:
        features = ad.obsm[use_rep]
        
    # Add nodes
    G.add_nodes_from((ad.obs.index.astype(int).tolist()), test=False, val=False, label=0)
    nx.set_node_attributes(G,dict(zip((ad.obs.index.astype(int).tolist()), features)), 'feature')
    # Add edges
    G.add_edges_from(res)
    nx.set_edge_attributes(G, dict(zip(res, d_res)), 'distance')
    
    # Remove components with less than N nodes
    N=args.min_conn_comp_size
    for component in tqdm(list(nx.connected_components(G))):
        if len(component)<N:
            for node in component:
                node_removed.append(node)
                G.remove_node(node)
    
    return G, ad[~ad.obs.index.astype(int).isin(node_removed),:].copy()


# One-hot-encoding
gene_list = df.gene.unique()
one_hot_encoding = dict(zip(gene_list,to_categorical(np.arange(gene_list.shape[0]),num_classes=gene_list.shape[0]).tolist()))
features = df['gene'].map(one_hot_encoding).tolist()
sample_ids = df.sample_id.unique()

# Create anndata object
adata = AnnData(np.array(features), obs=df)
adata.obsm["spatial"] = df.loc[:,['x','y','z']].values
adata.var_names = gene_list

# Compute spatial graph
if args.neigh_perc is not None:
    radius = computeNeighborDistance(adata, args.neigh_perc, spatial_key=args.spatial_key)
else:
    radius = args.radius

if len(sample_ids)>1:
    ad_l = []
    for s in adata.obs.sample_id.unique():
        adata_tmp = adata[adata.obs.sample_id==s].copy()
        sq.gr.spatial_neighbors( adata_tmp,
                         spatial_key=args.spatial_key, 
                         coord_type='generic',
                         radius=radius,
                         delaunay=args.delaunay,
                         key_added=args.spatial_key,
                         copy=False )
        
        ad_l.append(adata_tmp)
    adata = ann.concat(ad_l, uns_merge='same', pairwise=True)
else:
    sq.gr.spatial_neighbors( adata,
                         spatial_key=args.spatial_key, 
                         coord_type='generic',
                         radius=radius,
                         delaunay=args.delaunay,
                         key_added=args.spatial_key,
                         copy=False )
# Build stellargraph object
G, adata = buildGraph(adata, use_rep='X')
print("Average node degree: %f" + str(np.sum(list(dict(G.degree()).values()))/G.number_of_nodes()))

G = StellarGraph.from_networkx(
    G, edge_weight_attr="distance", node_features="feature"
)
print(G.info())

nodes = list(G.nodes())
number_of_walks = args.num_of_walks
length = args.walk_length
batch_size = args.batch_size
epochs = args.num_epochs
num_samples = [args.k1, args.k2]

unsupervised_samples = UnsupervisedSampler(G, nodes=nodes, length=length, number_of_walks=number_of_walks, seed=args.seed)

generator = GraphSAGELinkGenerator(G, batch_size, num_samples, seed=args.seed, weighted=True)
train_gen = generator.flow(unsupervised_samples)

layer_sizes = [args.l1, args.l2]

graphsage = GraphSAGE(
    layer_sizes=layer_sizes, generator=generator, aggregator=AttentionalAggregator, bias=True, dropout=0.0, normalize="l2", kernel_regularizer='l1'
)

# Build the model and expose input and output sockets of graphsage, for node pair inputs:
x_inp, x_out = graphsage.in_out_tensors()

prediction = link_classification(
    output_dim=1, output_act="sigmoid", edge_embedding_method='ip'
)(x_out)

earlystop_callback = tf.keras.callbacks.EarlyStopping(monitor='loss', mode='min', verbose=1, patience=1)

model = keras.Model(inputs=x_inp, outputs=prediction)
model.compile(
    optimizer=keras.optimizers.Adam(lr=0.5e-4),
    loss=keras.losses.binary_crossentropy,
    metrics=[keras.metrics.binary_accuracy]
)
model.summary()

multiprocessing=False if args.n_cpus==1 else True

# Start training
history = model.fit(
    train_gen,
    epochs=epochs,
    verbose=1,
    use_multiprocessing=multiprocessing,
    workers=args.n_cpus,
    shuffle=True,
    callbacks=[earlystop_callback]
)

# Extracting node embeddings
x_inp_src = x_inp[0::2]
x_out_src = x_out[0]
embedding_model = keras.Model(inputs=x_inp_src, outputs=x_out_src)

#Save the model
embedding_model.save("{}".format(FILE_PATH_MODEL_BASENAME))

nodes = list(G.nodes())
node_gen = GraphSAGENodeGenerator(G, batch_size, num_samples, seed=args.seed).flow(nodes)

node_embeddings = embedding_model.predict(node_gen, workers=args.n_cpus, verbose=1, use_multiprocessing=multiprocessing)
adata.obsm["spage2vec"] = node_embeddings
################################################################################

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

