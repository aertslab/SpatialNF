## Example data sets

#### `mouse_cortex_molecular_cartography`

This folder contains Molecular Cartography spatial transcriptomics data of 100 genes for three mouse cortex samples. Cells were segmented based on DAPI images. Counts of each transcript were computed for each segmented cell. Spatial data is provided in CSV compatible with SpatialNF. This data set was used in following publication:

_Bravo González-Blas, De Winter, Hulselmans, Hecker, Matetovici, Christiaens, ... & Aerts (2023). SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. Nat. Methods, 20:1355–1367._


#### `reference_scRNAseq`

This folder contains reference scRNAseq data that can be used for label-transfer to annotate spatial transcriptomics data:

- `Kleshchevnikov2022_mouseCortex_scRNAseq_98genes.h5ad` contains mouse brain scRNAseq data subset to cortical cell types. To reduce the size of the file, it was subset to genes that are shared with the `mouse_cortex_molecular_cartograph` data set. The data was obtained from following publication:
_Kleshchevnikov, V., Shmatko, A., Dann, E., Aivazidis, A., King, H. W., Li, T., ... & Bayraktar, O. A. (2022). Cell2location maps fine-grained cell types in spatial transcriptomics. Nat. Biotechnol., 40:661-671._
