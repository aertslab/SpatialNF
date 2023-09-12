# SpatialNF

## Spatial transcriptomics NextFlow pipelines

SpatialNF is a collection of Nextflow DSL2 pipelines for analyzing spatial transcriptomics data.

We offer pipelines for:
- basic processing of spatial transcriptomics data
- identification of spatially variable genes
- segmentation-free analysis
- label-transfer from single cell RNAseq to spatial transcriptomics data

SpatialNF is implemented in the VSN-framework: https://github.com/vib-singlecell-nf/vsn-pipelines.

### Supported data

SpatialNF can be used for FISH-based data like MERFISH or Molecular cartography and sequencing-based data like 10X Visium. As we do not offer automated segmentation pipelines within SpatialNF, FISH-based data has to be segmented in advance and converted into a supported data format.

For all data types, input data should contain raw counts. SpatialNF supports the following data formats:

| Data type  | Description |
| ------------- | ------------- |
| AnnData `.h5ad`  | AnnData object should contain an `.obsm` entry 'X_spatial' or 'spatial' storing coordinates of segmented cells or spots. Currently, our Docker images only support anndata <= 0.78. |
| 10X Spaceranger output | `outs` folder should contain the default 10X Spaceranger output |
| Spatial CSV files | A folder containing a coordinate file`coords.csv` and a count matrix file `matrix.csv`|
| Coordinates CSV files | A CSV file contaning coordinates of each transcript per row |

#### Spatial CSV file formats

A coordinate CSV file `coords.csv` contains three columns: an ID for the `parent_cell`, `X` and  `Y` coordinates:

```
parent_cell,X,Y
1,1952.8673508171798,3899.5127328012163
2,1946.5086419753086,4047.905679012346
3,1952.432242022379,3966.5445503522587
4,1963.7581227436824,4089.649097472924
5,1988.4492753623188,4047.0072463768115
...
```

A count matrix CSV file `matrix.csv` contains a column with `parent_cell` IDs matching IDs in the `coords.csv` and columns for each transcript:

```
parent_cell,Act79B,Act88F,AkhR,AstC-R2,Awh,CCAP-R,CG32121,Cralbp,FASN2
1,0,0,1,0,0,0,0,0,0
2,0,0,0,0,0,0,0,0,31
3,0,0,1,0,0,0,0,0,34
4,0,0,0,0,0,0,0,0,7
5,0,0,0,0,0,0,0,0,1
7,0,0,0,0,0,0,0,0,4
8,0,0,0,0,0,0,0,0,1
9,0,0,2,0,0,0,0,0,0
11,0,0,0,0,0,0,0,1,0
...
```

#### Coordinates CSV file formats

Coordinates CSV files contain the coordinates of each detected transcript per row. The should include a `header`, `x` and `y` columns.

```
gene,x,y
Rora,755,935
Rora,829,574
Rora,1071,1941
...
Slc17a7,2110,1458
Slc17a7,2110,1873
Slc17a7,2111,302
```
Counts per gene will be collated in a grid. The bin size can be speficied with `binsize`.


### Pipelines for basic processing

These pipelines assemble spatial transcriptomics data into AnnData and SCope (https://scope.aertslab.org) compatible loom files. They perform QC, filtering and clustering of the data

| Pipeline / entry point  | Description |
| ------------- | ------------- |
| single_sample | Process samples seperately |
| multi_sample | Compile and process samples together |

### Pipelines for identifying spatially variable genes

For detecting spatially variable genes, we implemented a pipeline using SpatialDE inlcuding their AEH approach for identifying spatial patterns.

| Pipeline / entry point  | Description |
| ------------- | ------------- |
| [single_sample_spatialde](https://github.com/aertslab/SpatialNF/tree/main/examples/spatialde) | Run `single_sample` pipeline and identify spatially variable genes and spatial patterns |
| multi_sample_spatialde | Run `multiple_sample` pipeline and identify spatially variable genes and spatial patterns |
| spatialde | Only run SpatialDE pipeline; input should be an AnnData object created by SpatialNF |


### Pipelines for label-transfer from scRNAseq

For label-transfer from scRNA-seq to spatial transcriptomics data, we offer pipelines for spot-based or segmented data using Tangram and SpaGE.
In addition, for segmentation-free label-transfer, SpatialNF contains a spage2vec pipeline.
Optionally, Squidpy can be used for computing enrichments of co-localized labels.
And Tangram can also be used to project gene expression from single cell data which will overwrite the count matrix.

The reference scRNAseq data should be a processed and filtered AnnData object `.h5ad` and contain raw counts, as well as the annotation as `obs` entry.

| Pipeline / entry point  | Description |
| ------------- | ------------- |
| [single_sample_tangram](https://github.com/aertslab/SpatialNF/tree/main/examples/tangram) | Run `single_sample` pipeline and Tangram for label-transfer. |
| multi_sample_tangram | Run `multiple_sample` pipeline and and Tangram for label-transfer. |
| tangram | Only run Tangram pipeline; input should be an AnnData object created by SpatialNF |
| single_sample_spage | Run `single_sample` pipeline and SpaGE for label-transfer. |
| multi_sample_spage | Run `multiple_sample` pipeline and and SpaGE for label-transfer. |
| spage2vec_spage_label_transfer | Perform neighborhood embedding analysis and label transfer with spage2vec. |


### Running SpatialNF pipelines

Initial configs can be generated with a `nextflow config` command. See the `SpatialNF/examples` folder for pipeline configuration files.
For example, to generate a config for mouse 10X spaceranger data and the `single_sample_spatialde` workflow:

```
nextflow config SpatialNF/main.nf \
    -profile mm10,tenx,singularity,single_sample,spatialde > single_sample_spatialde_aeh.config
```

After changing file names and parameters in the config, the pipeline can be run with the following command:

```
nextflow -C single_sample_spatialde_aeh.config run SpatialNF/main.nf \
   -entry single_sample_spatialde \
   -with-report report.html \
   -with-trace \
   -resume
```


### Output data

SpatialNF generates AnnData `.h5ad` and SCope (https://scope.aertslab.org) compatible loom `.loom` files.
Output data are written to `out/data/` including intermediate data files `out/data/intermediate/`.
Reports are stored as Jupyter notebooks and HTML files in `out/notebooks/`.


### Resource management

SpatialNF pipelines can be run locally or each step can be seperately submitted to as a job to a HPC.
Resource limits and parameters can be specified in the `process` section of a config file.


### Prerequisites

SpatialNF requires [singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html) to run Docker containers and [Nextflow](https://www.nextflow.io/). Currently, only Nextflow version 21.04 is supported. A compatible Netxtflow binary can be downloaded here: https://github.com/nextflow-io/nextflow/releases/download/v21.04.0/nextflow-21.04.0-all 


### Notes on Docker images
We are currenty providing Docker images at Docker hub on a free license. In case these Docker images become unavailable in the future, they can be rebuild from the `Dockerfile` in the workflow specific subfolders in the `src` directory.


## References

All pipelines are using SCANPY:

_Wolf, Angerer, & Theis (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biol., 19:1-5._
https://github.com/scverse/scanpy

SpatialDE:

_Svensson Teichmann & Stegle (2018). SpatialDE: identification of spatially variable genes. Nat. Methods, 15:343-346)_
https://github.com/Teichlab/SpatialDE

Tangram:

_Biancalani, Scalia, Buffoni, Avasthi, Lu, Sanger, ... & Regev (2021). Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat. Methods, 18:1352-1362._
https://github.com/broadinstitute/Tangram

SpaGE:

_Abdelaal, Mourragui, Mahfouz, & Reinders (2020). SpaGE: spatial gene enhancement using scRNA-seq. Nucleic Acids Res., 48:e107-e107._
https://github.com/tabdelaal/SpaGE

spage2vec:

_Partel & Waehlby (2021). Spage2vec: Unsupervised representation of localized spatial gene expression signatures. FEBS J., 288:1859-1870._
https://github.com/wahlby-lab/spage2vec

Squidpy:

_Palla, Spitzer, Klein, Fischer, Schaar, Kuemmerle, ... & Theis (2022). Squidpy: a scalable framework for spatial omics analysis. Nat. Methods, 19:171-178._
https://github.com/scverse/squidpy
