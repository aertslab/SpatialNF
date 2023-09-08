#!/usr/bin/env python3

import argparse
import os
import re
import scanpy as sc
import pandas as pd
import numpy as np
import loompy as lp
import anndata as ann
from scipy.sparse import csr_matrix
import json

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


in_formats = [
    '10x_cellranger_mex',
    '10x_cellranger_h5',
    '10x_spaceranger_visium',
    'spatial_csv'
    'tsv',
    'csv'
]

out_formats = [
    'h5ad'
]

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=str,
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "-i", "--input-format",
    action="store",
    dest="input_format",
    default="",
    help="Input format of the file to be converted. Choose one of: {}.".format(', '.join(in_formats))
)

parser.add_argument(
    "-s", "--sample-id",
    type=str,
    dest="sample_id",
    default=None,
    action='store'
)

parser.add_argument(
    "-k", "--group-name",
    type=str,
    dest="group_name",
    default=None,
    action='store',
    help="Name of the group which the given input files are from. A new column named by this group_name will be added in anndata.obs"
)

parser.add_argument(
    "-l", "--group-value",
    type=str,
    dest="group_value",
    default=None,
    action='store',
    help="Value of the group which the given input files are from. The group_name column in anndata.obs will be populated with this group_value."
)

parser.add_argument(
    "-t", "--tag-cell-with-sample-id",
    type=str2bool,
    action="store",
    dest="tag_cell_with_sample_id",
    default=False,
    help="Tag each cell with the given sample_id."
)

parser.add_argument(
    "-r", "--remove-10x-gem-well",
    type=str2bool,
    action="store",
    dest="remove_10x_gem_well",
    default=False,
    help="If tag_cell_with_sample_id is passed, remove the GEM well number from the barcode."
)

parser.add_argument(
    "-u", "--make-var-index-unique",
    type=str2bool,
    action="store",
    dest="make_var_index_unique",
    default=False,
    help="Make the var index unique of the AnnData."
)

parser.add_argument(
    "-w", "--use-raw",
    type=str2bool,
    action="store",
    dest="use_raw",
    default=False,
    help="Replace X entry with raw.X matrix. This only is used when input format is h5ad."
)

parser.add_argument(
    "-o", "--output-format",
    action="store",  # optional because action defaults to "store"
    dest="output_format",
    default="",
    help="Output format which the file should be converted to. Choose one of: {}.".format(', '.join(out_formats))
)

parser.add_argument(
    "--scale-spatial",
    action="store",  # optional because action defaults to "store"
    dest="spatial_scale",
    type=float,
    default=20,
    help="Real number to scale 'X_spatial' to for SCope. Set to 0 for no scaling."
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]
INPUT_FORMAT = args.input_format
OUTPUT_FORMAT = args.output_format

def check_spatial_csv_path(path):
    # Sanity checks
    if not os.path.isdir(path):
        raise Exception(
            "Expecting a directory with coords.csv and matrix.csv files when converting from spatial_csv. {} does not seem to be one.".format(
                path
            )
        )
    if not os.path.exists(path):
        raise Exception("VSN ERROR: The given directory {} does not exist.".format(path))
    if (
        not os.path.exists(os.path.join(path, "coords.csv")) or not os.path.exists(os.path.join(path, "matrix.csv"))
    ):
        raise Exception(
            "Expecting a directory with coords.csv and matrix.csv files when converting from spatial_csv. The given directory {} is not a proper folder. No csv file[s] found.".format(
                path
            )
        )


def check_10x_cellranger_mex_path(path):
    # Sanity checks
    if not os.path.isdir(path):
        raise Exception(
            "Expecting a directory with an .mtx file when converting from 10x_cellranger_mex. {} does not seem to be one.".format(
                path
            )
        )
    if not os.path.exists(path):
        raise Exception("VSN ERROR: The given directory {} does not exist.".format(path))
    if not (
        not os.path.exists(os.path.join(path, "matrix.mtx")) or not os.path.exists(os.path.join(path, "matrix.mtx.gz"))
    ):
        raise Exception(
            "The given directory {} is not a proper 10xGenomics CellRanger folder. No .mtx[.gz] file found.".format(
                path
            )
        )

    
def check_10x_spaceranger_visium_path(path):
    # Sanity checks

    # expected image files in 'spatial' sub dir
    image_files = ["tissue_hires_image.png", "tissue_lowres_image.png", "aligned_fiducials.jpg", "scalefactors_json.json", "detected_tissue_image.jpg", "tissue_positions_list.csv"]
    
    if not os.path.isdir(path):
        raise Exception(
            "Expecting a directory with an 'filtered_feature_bc_matrix.h5' file and 'spatial' sub dir containing image data when converting from 10x_spaceranger_visium. {} does not seem to be one.".format(
                path
            )
        )
    if not os.path.exists(path):
        raise Exception("VSN ERROR: The given directory {} does not exist.".format(path))

    # check filtered feature bit counts file
    if not os.path.exists(os.path.join(path, "filtered_feature_bc_matrix.h5")):
        raise Exception(
            "The given directory {} is not a proper 10x Genomics Space Ranger folder. No 'filtered_feature_bc_matrix.h5' file found.".format(
                path
            )
        )
    
    # check spatail subdir
    image_path = os.path.join(path, "spatial")
    
    if not os.path.isdir(image_path):
        raise Exception(
            "Expecting a 'spatial'sub dir containing image data when converting from 10x_spaceranger_visium. {} does not seem to be one.".format(
                image_path
            )
        )
    if not os.path.exists(image_path):
        raise Exception("VSN ERROR: The image directory {} does not exist.".format(image_path))
    # check image files 
    for img_file in image_files:
        if not os.path.exists(os.path.join(image_path, img_file)):
            raise Exception(
                "The given directory {} is not a proper 10x Genomics Space Ranger image folder. No {} file found.".format(
                    image_path, img_file
                )
            )

    
    

def add_obs_column(adata, column_name, value):
    # Annotate the file with the sample ID
    print(f"Adding new obs column '{column_name}' with value '{value}'...")
    adata.obs[column_name] = value
    return adata


def tag_cell(adata, tag, remove_10x_gem_well=False):
    # Check the number of untagged cells
    # We consider an untagged cell as matching the following pattern: [barcode-id]-[sample-index] where
    # - [barcode-id] is sequence of A,C,G,T letters
    # - [sample-index] is a natural number
    if remove_10x_gem_well:
        num_untagged_cells = sum(list(map(lambda x: len(re.findall(r"[ACGT]*-[0-9]+$", x)), adata.obs.index)))
        if num_untagged_cells != 0:
            print("Appending sample ID to cell barcode with 10x GEM well removal...")
            adata.obs.index = list(map(lambda x: re.sub(r"([ACGT]*)-.*", rf'\1-{tag}', x), adata.obs.index))
    else:
        print("Appending sample ID to cell barcode...")
        adata.obs.index = [cell_barcode + "___" + tag for cell_barcode in adata.obs.index]
    return adata


def update_obs(adata, args):
    # Add sample ID
    adata = add_obs_column(
        adata=adata,
        column_name="sample_id",
        value=args.sample_id
    )
    # If tag_cell_with_sample_id is given, add the sample ID as suffix
    if args.tag_cell_with_sample_id:
        adata = tag_cell(
            adata=adata,
            tag=args.sample_id,
            remove_10x_gem_well=args.remove_10x_gem_well
        )
    # Add group_value as obs entry with group_name as column name
    if args.group_name is not None and args.group_value is not None:
        adata = add_obs_column(
            adata=adata,
            column_name=args.group_name,
            value=args.group_value
        )
    return adata


def update_var(adata, args):
    adata.var.index = adata.var.index.astype(str)
    # Check if var index is unique
    if len(np.unique(adata.var.index)) < len(adata.var.index) and not args.make_var_index_unique:
        raise Exception("VSN ERROR: AnnData var index is not unique. This can be fixed by making it unique. To do so update the following param 'makeVarIndexUnique = true' (under params.utils.sc_file_converter) in your config.")
    if len(np.unique(adata.var.index)) < len(adata.var.index) and args.make_var_index_unique:
        adata.var_names_make_unique()
        print("Making AnnData var index unique...")
    return adata


if INPUT_FORMAT == '10x_cellranger_mex' and OUTPUT_FORMAT == 'h5ad':
    check_10x_cellranger_mex_path(path=FILE_PATH_IN)
    # Convert
    print("Reading 10x data from MEX format...")
    adata = sc.read_10x_mtx(
        FILE_PATH_IN,  # the directory with the `.mtx` file
        var_names='gene_symbols',  # use gene symbols for the variable names (variables-axis index)
        cache=False
    )
    # Add/update additional information to observations (obs)
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == '10x_cellranger_h5' and OUTPUT_FORMAT == 'h5ad':
    if not os.path.exists(FILE_PATH_IN):
        raise Exception("VSN ERROR: The given file {} does not exist.".format(FILE_PATH_IN))
    # Convert
    print("Reading 10x data from HDF5 format...")
    adata = sc.read_10x_h5(
        FILE_PATH_IN
    )
    # Add additional information
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT in ['tsv', 'csv'] and OUTPUT_FORMAT == 'h5ad':
    if INPUT_FORMAT == 'tsv':
        delim = '\t'
    elif INPUT_FORMAT == 'csv':
        delim = ','
    # Expects csv/tsv to have features as rows and observations as columns
    adata = sc.read_csv(
        FILE_PATH_IN,
        delimiter=delim,
        first_column_names=True
    ).T
    # Convert to sparse matrix
    adata.X = csr_matrix(adata.X)
    # Add additional information
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == 'spatial_csv' and OUTPUT_FORMAT == 'h5ad':
    check_spatial_csv_path(path=FILE_PATH_IN)
    # Convert
    print("Reading spatial csv format...")
    # Coordinates - Expects csv to have coordinate axes as columns and observations as rows
    coords_df = pd.read_csv(
        os.path.join(FILE_PATH_IN, "coords.csv"),
        sep=',',
        index_col=0
    )
    # Data matrix - Expects csv to have features as columns and observations as rows
    adata = sc.read_csv(
        os.path.join(FILE_PATH_IN, "matrix.csv"),
        delimiter=',',
        first_column_names=True
    )

    # assign 'Gene' and 'CellID' 
    adata.var["Gene"] = adata.var_names
    adata.obs["CellID"] = adata.obs_names

    # assign spatial info
    adata.obsm['spatial'] = coords_df.values
    # convert spatial coords into float
    adata.obsm['X_spatial'] = np.float32(adata.obsm['spatial'])[:,:2]
    # center data for SCope
    avg_coords = adata.obsm['X_spatial'].sum(axis=0)/adata.obsm['X_spatial'].shape[0]
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] - avg_coords
    # scale data for SCope, default of 20 such that X axis is between [-10,10]
    if args.spatial_scale > 0:
        scfactor = args.spatial_scale / (np.max(adata.obsm['X_spatial'][:,0])-np.min(adata.obsm['X_spatial'][:,0]))
    else:
        scfactor = 1
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] * scfactor

    # Convert to sparse matrix
    adata.X = csr_matrix(adata.X)
    # Add additional information
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]

    print("Writing spatial data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == 'h5ad' and OUTPUT_FORMAT == 'h5ad':
    adata = sc.read_h5ad(
        FILE_PATH_IN
    )

    if args.use_raw and adata.raw is not None:
        # Replace X with raw.X
        print("Replacing X with raw.X...")
        adata = sc.AnnData(
            X=adata.raw.X,
            obs=adata.obs,
            var=adata.raw.var
        )

    # Add additional information
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == '10x_spaceranger_visium' and OUTPUT_FORMAT == 'h5ad':
    check_10x_spaceranger_visium_path(path=FILE_PATH_IN)
    # Convert
    print("Reading 10x Visium data...")
    adata = sc.read_visium(
        FILE_PATH_IN,
        count_file='filtered_feature_bc_matrix.h5',
        load_images=True)

    # assign 'Gene' and 'CellID' 
    adata.var["Gene"] = adata.var_names
    adata.obs["CellID"] = adata.obs_names

    # convert spatial info into float
    adata.obsm['X_spatial'] = np.float32(adata.obsm['spatial'])
    
    # Add/update additional information to observations (obs)
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    adata = adata[:, np.sort(adata.var.index)]

    # center data for SCope
    avg_coords = adata.obsm['X_spatial'].sum(axis=0)/adata.obsm['X_spatial'].shape[0]
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] - avg_coords
    # scale data for SCope, default of 20 such that X axis is between [-10,10]
    if args.spatial_scale > 0:
        scfactor = args.spatial_scale / (np.max(adata.obsm['X_spatial'][:,0])-np.min(adata.obsm['X_spatial'][:,0]))
    else:
        scfactor = 1
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] * scfactor
    
    print("Writing 10x data to h5ad...")
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))

elif INPUT_FORMAT == 'loom' and OUTPUT_FORMAT == 'h5ad':

    # process loom
    loom = lp.connect(FILE_PATH_IN, validate=False)

    ex_mtx = pd.DataFrame(loom[:, :], index=loom.ra.Gene, columns=loom.ca.CellID).T
    col_attrs = {k: v for k, v in loom.ca.items()}
    row_attrs = {k: v for k, v in loom.ra.items()}
    global_attrs = {k: v for k, v in loom.attrs.items()}
    dict_metadata = json.loads(loom.__dict__['attrs']['MetaData'])

    # close loom
    loom.close()

    # create anndata
    adata = ann.AnnData(X=np.float32(ex_mtx.to_numpy()))
    adata.obs_names = ex_mtx.index
    adata.var_names = ex_mtx.columns

    # add column attributes
    for key, value in row_attrs.items():
        if not isinstance(value[0], np.void):
            adata.var[key] = value        
    # add column attributes
    for key, value in col_attrs.items():
        if not isinstance(value[0], np.void):
            adata.obs[key] = value
    # add column attributes
    for key, value in global_attrs.items():
        adata.uns[key] = value
    
    # add clustering
    if 'clusterings' in dict_metadata:

        for clusidx in range(len(dict_metadata['clusterings'])):

            clustering_name = dict_metadata['clusterings'][clusidx]['name']

            # get cluster id --> name mapping; use first annotation if present otherwise description
            map_clus_names = {}
            for clus in dict_metadata['clusterings'][clusidx]['clusters']:
                if 'cell_type_annotation' in clus:
                    map_clus_names[clus['id']] = clus['cell_type_annotation'][0]['data']['annotation_label'].replace('/', '_')
                else:
                    map_clus_names[clus['id']] = clus['description'].replace('/', '_')

            # add first clustering with ranked genes group, others only as obs entry
            if clusidx == 0:
                clustering_algorithm = dict_metadata['clusterings'][0]['group']
            
                match_res = re.search(r'[0-9]+\.?[0-9]?$', clustering_name)
                if match_res:
                    clustering_resolution = re.search(r'[0-9]+\.?[0-9]?$', clustering_name)[0]
                else:
                    clustering_resolution = 0
                mmethod = re.search(r'Average log fold change from (.*)', dict_metadata['clusterings'][0]['clusterMarkerMetrics'][0]['description'])
            
                if mmethod:
                    cluster_marker_method = mmethod[1]
                else:
                    cluster_marker_method = "None"
            
                if cluster_marker_method != "None":
                    adata.uns["rank_genes_groups"] = {}
                    adata.uns["rank_genes_groups"]["params"] = {}
                    adata.uns["rank_genes_groups"]["params"]["groupby"] = clustering_algorithm
                    adata.uns[clustering_algorithm] = {}
                    adata.uns[clustering_algorithm]["params"] = {}
                    adata.uns[clustering_algorithm]["params"]["resolution"] = clustering_resolution
                    adata.uns["rank_genes_groups"]["params"]["method"] = cluster_marker_method
                
                    # add marker genes
                    # init empty dict
                    adata.uns["rank_genes_groups"]['names'] = {}
                    adata.uns["rank_genes_groups"]['pvals_adj'] = {}
                    adata.uns["rank_genes_groups"]['logfoldchanges'] = {}
                
                    for clusid in map_clus_names.values():
                        adata.uns["rank_genes_groups"]['names'][clusid] = []
                        adata.uns["rank_genes_groups"]['pvals_adj'][clusid] = []
                        adata.uns["rank_genes_groups"]['logfoldchanges'][clusid] = []
                
                    # get marker genes
                    for i in range(len(row_attrs['ClusterMarkers_0'])):
                        gene = row_attrs['Gene'][i]
                        for idx, e in enumerate(row_attrs['ClusterMarkers_0'][i]):
                            if e != 0:
                                pval = row_attrs['ClusterMarkers_0_pval'][i][idx]
                                logfc = row_attrs['ClusterMarkers_0_avg_logFC'][i][idx]
                                clusid = map_clus_names[idx]
                                adata.uns["rank_genes_groups"]['names'][clusid].append(gene)
                                adata.uns["rank_genes_groups"]['pvals_adj'][clusid].append(pval)
                                adata.uns["rank_genes_groups"]['logfoldchanges'][clusid].append(logfc)
                
                # add obs entry for clustering
                adata.obs[clustering_algorithm] = [ map_clus_names[clus[clusidx]] if clus[clusidx] in map_clus_names.keys() else clus[clusidx] for clus in col_attrs['Clusterings'] ]

            adata.obs[clustering_algorithm] = [ map_clus_names[clus[clusidx]] if clus[clusidx] in map_clus_names.keys() else clus[clusidx] for clus in col_attrs['Clusterings'] ]
    
    # add embeddings
    map_embedding_names = {
        'HVG UMAP': 'X_umap',
        'HVG t-SNE': 'X_tsne',
        'HVG PC1/PC2': 'X_pca',
        'Spatial': 'X_spatial'
        }

    if 'embeddings' in dict_metadata:

        # initialize embeddings        
        dict_embed_x = {}
        dict_embed_y = {}
        n_spots = len(col_attrs['Embeddings_X'])
        veczero = np.zeros(n_spots)
        embeds = [ embed['name'] for embed in dict_metadata['embeddings'] ]
        n_embeds = len(embeds)
        
        for embed in embeds:
            dict_embed_x[embed] = veczero.copy()
            dict_embed_y[embed] = veczero.copy()
            
        for i in range(n_spots):
            for j in range(n_embeds):
                dict_embed_x[embeds[j]][i] = col_attrs['Embeddings_X'][i][j]
                dict_embed_y[embeds[j]][i] = col_attrs['Embeddings_Y'][i][j]
        
        for embed in embeds:
            if embed in map_embedding_names:
                name_embed = map_embedding_names[embed]
            else:
                name_embed = "X_" + embed
            
            x = dict_embed_x[embed]
            y = dict_embed_y[embed]
            
            adata.obsm[name_embed] = np.float32(np.vstack((x, y)).T)
    
    
    # Add additional information
    adata = update_obs(adata=adata, args=args)
    # Add/update additional information to features (var)
    adata = update_var(adata=adata, args=args)
    # Sort var index
    #adata = adata[:, np.sort(adata.var.index)]
    adata.write_h5ad(filename="{}.h5ad".format(FILE_PATH_OUT_BASENAME))
else:
    raise Exception(
        "File format conversion {0} --> {1} hasn't been implemented yet.".format(INPUT_FORMAT, OUTPUT_FORMAT))

print("Done!")
