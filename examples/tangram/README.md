## Tangram workflows

### Single sample Tangram

As an example, we provide a config file `single_sample_tangram.config` for Molecular Cartography spatial transcriptomics data: `SpatialNF/examples/data/mouse_cortex_molecular_cartography`

A config file for a  single sample Tangram workflow for input as CSV files can be generated with following command:

```
nextflow config SpatialNF/main.nf -profile profile mm10,singularity,spatial_csv,single_sample,tangram,squidpy > single_sample_tangram.config
```
where `SpatialNF` is the path to the SpatialNF parent directory. This config file should be adjusted to contain the path to the data and proper parameters.


#### input data

The input data is specified in the `data` section. For label transfer, a reference AnnData object containing the annotation and the spatial transcriptomics data should be specified. The reference AnnData object can be found at `SpatialNF//data/reference_scRNAseq/Kleshchevnikov2022_mouseCortex_scRNAseq_98genes.h5ad`. In this case, we use CSV files as spatial transcriptomics input data:

```
   data {
      spatial_csv = '/lustre1/project/stg_00002/lcb/nhecker/software/git/SpatialNF/examples/data/mouse_cortex_molecular_cartography/A1/outs'
      reference_data {
         file_path = '/lustre1/project/stg_00002/lcb/nhecker/software/git/SpatialNF/examples/data/reference_scRNAseq/Kleshchevnikov2022_mouseCortex_scRNAseq_98genes.h5ad'
         suffix = 'h5ad'
      }
   }
```
Always use full file paths in order to allow singularity to find them. The reference AnnData object should contain filtered data, raw counts, and the annation as `.obs` entry.

The `data` parameters are specific to the input data type. For 10X spaceranger output as input data, use a `spaceranger_visium` entry:

```
   data {
       tenx {
          spaceranger_visium = '/lustre1/data/spatial/Visium_brain/spaceranger/*/outs'
       }
       ...
   }
```
where `*` is a wildcard to use all subfolders; one sample per subfolder. If the input data are already AnnData objects, following data entry can be used:

```
   data {
        h5ad {
         file_paths = '/lustre1/data/spatial/mouse_brain_*.h5ad'
         suffix = 'h5ad' 
        }
	...
   }
```

If the input data is given as individual transcripts, they will be binned in a grid. The bin size can be specified with `binsize`:

```
   data {
      coordinates_csv {
         file_paths = '/lustre1/project/stg_00002/lcb/nhecker/software/git/SpatialNF/examples/data/mouse_cortex_molecular_cartography/*/outs/*_coordinates.csv'
         suffix = 'csv'
      }
      binsize = 200
      ...
   }
```


Singularity parameters should be adjusted to contain the mount points where the data is stored:

```
singularity {
   enabled = true
   autoMounts = true
   runOptions = '--cleanenv -H $PWD -B ${HOME},/lustre1,/scratch,/staging'
}
```

#### Filtering

For smFISH based data, filtering criteria should be changed to less strict values. A minimum of 3 genes per cell and 3 cells per gene is a reasonable value for a data set of 100 genes:

```
   tools {
      scanpy {
      ...
         filter {
            report_ipynb = '/src/scanpy/bin/reports/sc_filter_qc_report.ipynb'
            cellFilterStrategy = 'fixedthresholds'
            cellFilterMinNGenes = 3
            cellFilterMaxNGenes = 4000
            cellFilterMaxPercentMito = 0.15
            geneFilterMinNCells = 3
            off = 'h5ad'
            outdir = 'out'
         }
```

#### Tangram parameters:

Tangram parameters are specified in the `tangram` parameters section:

```
      tangram {
         container = 'nhecker/tangram-sc:1.0.2'
         off = 'h5ad'
	 // name of cell type annotation used for mapping
         annotation = 'BRAIN_cell_type'
	 // 'all' genes for MERFISH-like data, 'marker_genes', or 'list' for a user-specified gene list
         gene_selection_method = 'all'
	 // method for computing marker genes if not computed already
         rank_gene_method = 'wilcoxon'
         number_genes = 100
         qvalue = 0.05
	 // maximum sparsity per cluster
         max_sparsity = 0.5
         mapping_mode = 'cells'
	 // log-normalize scRNAseq data
         normalize_scRNAseq = true
	 // log-normalize spatial data
         normalize_spatial = true
	 // change to true to project gene expression from single cell reference data
         project_gex = false
	 // device for computing mapping, e.g. '0' for 'cuda: 0' change to 'cpu' for using CPU; 'any' will let torch pick GPU
         device = 'any'
         report_ipynb = '/src/tangram/bin/reports/tangram_celltype_projection_report.ipynb'
	 // normalize scores after mapping cell types
         normalize_celltype_scores = true
	 // mapping score quantile for binary assignments of celltypes
         quantile_mapping_score = 0.95
	 // compute squidpy neighborhood statistics, useful cell-segmented spatial data (set to 'true')
         squidpy_statistics = true
      }
```

`annotation` specifies the labels to be transferred from scRNAseq to spatial transcriptomics data.
Tangram computes a mapping between each single cell of the reference scRNAseq to spot/segmented cell. This provides a continuous score per transferred label for each cell. As often proportions between cells in single cell data and spatial transcriptomics data do not, match normalizing scores by the maximum score per cell type can be usefule `normalize_celltype_scores = true`.

For sequencing-based technologies like 10X Visium or Stereoseq, a subset of marker genes should be used instead of all shared genes between scRNAseq and spatial transcriptomics data. How genes are selected can be specified by `gene_selection_method`:
- `all`: all shared genes
- `list`: specify the full path to a file with a list of manually picked genes as `gene_list`: one gene per row. 
- `marker_genes`: marker genes computed per cell type are computed by SCANPY
-- `rank_gene_method`: specifies the test to be used for identifying marker genes
-- `number_genes`: specifies the N top marker genes per cell type to be used
-- `qvalue`: specifies the q-value threshold for genes to be selected
-- `max_sparsity`: specifies the sparsity threshold for a gene to be considered as marker gene

To project gene expression from scRNAseq to single cell data, set `project_gex = true`.

To assign a labels to spots/segmented cells the highest scoring label will be picked this annotation will be specified as `.obs` entry `tangram_best_` + annotation, e.g. `tangram_best_BRAIN_cell_type`. Often it might be more useful to assign labels based on a score quantile like the 95%-quantile. This quantile can be changed with `quantile_mapping_score` and will be added as seperate `.obs` entry.

Squidpy statistics or label neighnorhood enrichments can be computed. Sometimes this will cause issues. It can be switched by setting `squidpy_statistics = false`.


Additional options:
```
    tools {
        tangram {               
		...
                // downsample to maximum number of cells per cluster
                //max_cells = 1000
                // user-specified gene list file, one gene per row, no header
                //gene_list = ''
        }
    }
```

As Tangram computes a matrix N spots x M single cells, it can require large amounts of videomemory. To limit the amount of required videomemory, cells can be downsampled to a specific maximum number of cells per cell type with `max_cells = 1000`.


See https://github.com/broadinstitute/Tangram for more information on Tangram.


#### Squidpy statistics

Parameters for squidpy statistics are specified in the `squidpy` section:

```
      squidpy {
         container = 'gapartel/squidpy:1.1.2.dev'
         report {
            annotations_to_plot = ['tangram_best_BRAIN_cell_type']
         }
         spatial_graph {
            spatial_key = 'spatial'
            coord_type = 'grid'
            off = 'h5ad'
         }
         spatial_statistics {
            report_ipynb = '/src/squidpy/bin/reports/spatial_squidpy_statistics.ipynb'
            cluster_key = 'tangram_best_BRAIN_cell_type'
            connectivity_key = 'spatial'
            off = 'h5ad'
         }
         lr_analysis {
            cluster_key = 'tangram_best_BRAIN_cell_type'
            off = 'h5ad'
         }
      }
```

Entries that should be adjust are the `annotations_to_plot = ['tangram_best_BRAIN_cell_type']` and `cluster_key = 'tangram_best_BRAIN_cell_type'` entries.

See https://github.com/scverse/squidpy for more information on squidpy.


#### Running Tangram workflows

After adjusting the config file, following command can be used to execute the workflow:

```
nextflow -C single_sample_tangram.config run SpatialNF/main.nf -entry single_sample_tangram -resume -with-report
```

To compute mappings for all samples together, tangram can be run as multiple sample workflow:

```
nextflow -C multi_sample_tangram.config run SpatialNF/main.nf -entry multi_sample_tangram -resume -with-report
```

If AnnData objects for spatial transcriptomics data were already generated by SpatialNF, the single sample or multi sample workflow part can be skipped with:

```
nextflow -C tangram.config run SpatialNF/main.nf -entry tangram -resume -with-report
```

Input data should be already filtered AnnData in that case containing an `.obsm['X_spatial']` entry for the coordinates and raw counts.

```
   data {
       h5ad {
         file_paths = '/data/spatial/mouseCortex/A1_raw_filtered.h5ad'
         suffix = 'h5ad' 
       }
       reference_data {
         file_path = '/lustre1/project/stg_00002/lcb/nhecker/software/git/SpatialNF/examples/data/reference_scRNAseq/Kleshchevnikov2022_mouseCortex_scRNAseq_98genes.h5ad'
         suffix = 'h5ad'
      }
   }
```

#### Output

Tangram workflows generate:
- AnnData objects containing the projected labels e.g. `out/data/A1.TANGRAM_CELLTYPES.h5ad` and mapping `out/data/A1.TANGRAM_MAPPING.h5ad`
- SCope compatible Loom files e.g. `out/data/A1.TANGRAM_CELLTYPES_scope.loom`
- reports e.g. `out/notebooks/A1.TANGRAM_report.html`
