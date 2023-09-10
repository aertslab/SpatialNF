## SpatialDE workflows

### Single sample spatial DE

As an example, we provide a config file `single_sample_spatialde_aeh.config` for Molecular Cartography spatial transcriptomics data: `SpatialNF/examples/data/mouse_cortex_molecular_cartography`

A config file for a  single sample spatialDE workflow for input as CSV files can be generated with following command:

```
nextflow config SpatialNF/main.nf -profile mm10,singularity,spatial_csv,single_sample,spatialde > single_sample_spatialde_aeh.config
```
where `SpatialNF` is the path to the SpatialNF parent directory. This config file should be adjusted to contain the path to the data and proper parameters.


#### input data

The input data is specified in the `data` section. In this case, we use CSV files as input:

```
   data {
      spatial_csv = '/lustre1/project/stg_00002/lcb/nhecker/software/git/SpatialNF/examples/data/mouse_cortex_molecular_cartography/*/outs'
   }
```
where `*` is a wildcard to use all three subfolders; one subfolder per sample. Always use full file paths in order to allow singularity to find them.

The `data` parameters are specific to the input data type. For 10X spaceranger output as input data, use a `spaceranger_visium` entry:

```
   data {
       tenx {
          spaceranger_visium = '/lustre1/data/spatial/Visium_brain/spaceranger/*/outs'
       }
   }
```

If the input data are already AnnData objects, following data entry can be used:

```
   data {
        h5ad {
         file_paths = '/lustre1/data/spatial/mouse_brain_*.h5ad'
         suffix = 'h5ad' 
        }
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

#### SpatialDE parameters:

SpatialDE parameters are specified in the `spatialde` parameters section:

```
      spatialde {
         container = 'nhecker/spatialde:latest'
         off = 'h5ad'
         min_cells = 3
         thr_qval = 0.05
         method_pval_correction = 'bonferroni'
         min_fsv = 0.5
         normalize_naivede = true
         run_aeh = true
         c = 5
         estimate_l = 'mean'
         adjust_l = 0.2
         report_ipynb = '/src/spatialde/bin/reports/spatialde_spatially_variable_genes_report.ipynb'
         report_aeh_ipynb = '/src/spatialde/bin/reports/spatialde_automatic_expression_histology.ipynb'
      }
```
Filtering can be adjusted for SpatialDE `min_cells`. Significance for spatially variable genes is specified by q-value threshold `thr_qval` and a minimum fraction of spatial variance `min_fsv`.

SpatialDE uses Automatic Expression Histology (AEH) to group spatially variable genes into patterns. This can be very computationally demanding. To skip AEH, set `run_aeh = false`.
AEH is based on a length scale `l`. This length scale can be set manually or automatically selected by the mean value of of estimated length scales `estimate_l = 'mean'` plus a fraction of then mean value `adjust_l = 0.2`.
The number of patterns to find is defined by `c = 5`. 

For more details see: https://github.com/Teichlab/SpatialDE

#### Running SpatialDE workflow

After adjusting the config file, following command can be used to execute the workflow:

```
nextflow -C single_sample_spatialde_aeh.config run SpatialNF/main.nf -entry single_sample_spatialde -resume -with-report
```


#### Output

SpatialDE workflows generate:
- AnnData objects e.g. `out/data/A1.SPATIALDE__SPATIAL_VARIABLE_GENES.h5ad`; SpatialDE results added as `adata.uns['spatialDE']`, `adata.uns['spatialDE_AEH_histology_results']` and `adata.uns['spatialDE_AEH_patterns']`
- SCope compatible Loom files `out/data/A1.SPATIALDE_scope.loom`.
- reports about spatially variable genes e.g. `A1.spatialDE_variable_genes_report.html`
- reports about spatial pattern e.g. `A2.spatialDE_aeh_report.html `
