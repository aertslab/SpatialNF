nextflow.enable.dsl=2

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the utils module:
include {
    getBaseName;
} from '../src/utils/processes/files.nf'
include {
    clean;
    SC__FILE_CONVERTER;
    SC__FILE_CONCATENATOR;
} from '../src/utils/processes/utils.nf' params(params)
include {
    COMBINE_BY_PARAMS;
    PUBLISH;
} from '../src/utils/workflows/utils.nf' params(params)
include {
    FINALIZE;
} from '../src/utils/workflows/finalize.nf' params(params)
include {
    SC__H5AD_TO_FILTERED_LOOM;
} from '../src/utils/processes/h5adToLoom.nf' params(params)
include {
    FILE_CONVERTER;
} from '../src/utils/workflows/fileConverter.nf' params(params)
include {
    FILTER_AND_ANNOTATE_AND_CLEAN;
} from '../src/utils/workflows/filterAnnotateClean.nf' params(params)
include {
    UTILS__GENERATE_WORKFLOW_CONFIG_REPORT;
} from '../src/utils/processes/reports.nf' params(params)

////////////////////////////////////////////////////////
//  Import sub-workflows/processes from the tool module:
include {
    QC_FILTER;
} from '../src/scanpy/workflows/qc_filter.nf' params(params)
include {
    NORMALIZE_TRANSFORM;
} from '../src/scanpy/workflows/normalize_transform.nf' params(params)
include {
    HVG_SELECTION;
} from '../src/scanpy/workflows/hvg_selection.nf' params(params)
include {
    SC__SCANPY__REGRESS_OUT;
} from '../src/scanpy/processes/regress_out.nf' params(params)
include {
    NEIGHBORHOOD_GRAPH;
} from '../src/scanpy/workflows/neighborhood_graph.nf' params(params)
include {
    DIM_REDUCTION_PCA;
} from '../src/scanpy/workflows/dim_reduction_pca.nf' params(params)
include {
    DIM_REDUCTION_TSNE_UMAP;
} from '../src/scanpy/workflows/dim_reduction.nf' params(params)
include {
    SC__SCANPY__CLUSTERING_PARAMS;
} from '../src/scanpy/processes/cluster.nf' params(params)
include {
    CLUSTER_IDENTIFICATION;
} from '../src/scanpy/workflows/cluster_identification.nf' params(params)
include {
    SC__DIRECTS__SELECT_DEFAULT_CLUSTERING
} from '../src/directs/processes/selectDefaultClustering.nf'
// reporting:
include {
    SC__SCANPY__MERGE_REPORTS;
} from '../src/scanpy/processes/reports.nf' params(params)
include {
    SC__SCANPY__REPORT_TO_HTML;
} from '../src/scanpy/processes/reports.nf' params(params)


workflow multi_sample {

    take:
        data

    main:
        /*******************************************
        * Data processing
        */
        // To avoid Variable `params` already defined in the process scope
        def scanpyParams = params.tools.scanpy

        out = data | \
            SC__FILE_CONVERTER | \
            FILTER_AND_ANNOTATE_AND_CLEAN

        filtered = scanpyParams?.filter ? QC_FILTER( out ).filtered : out 

        if(params.utils?.file_concatenator) {
            concatenated = SC__FILE_CONCATENATOR( 
                filtered.map {
                    it -> it[1]
                }.toSortedList( 
                    { a, b -> getBaseName(a, "SC") <=> getBaseName(b, "SC") }
                ) 
            )
        } else {
            concatenated = filtered
        }

        transformed_normalized = scanpyParams?.data_transformation && scanpyParams?.normalization
            ? NORMALIZE_TRANSFORM( concatenated ) : concatenated

        out = HVG_SELECTION( transformed_normalized )
        DIM_REDUCTION_PCA( out.scaled )
        NEIGHBORHOOD_GRAPH( DIM_REDUCTION_PCA.out )
        DIM_REDUCTION_TSNE_UMAP( NEIGHBORHOOD_GRAPH.out )
        CLUSTER_IDENTIFICATION(
            NORMALIZE_TRANSFORM.out,
            DIM_REDUCTION_TSNE_UMAP.out.dimred_tsne_umap,
            "No Batch Effect Correction"
        )
        marker_genes = CLUSTER_IDENTIFICATION.out.marker_genes.map {
            it -> tuple(it[0], it[1], null)
        }

        // Finalize
        FINALIZE(
            params.utils?.file_concatenator ? SC__FILE_CONCATENATOR.out : SC__FILE_CONVERTER.out,
            CLUSTER_IDENTIFICATION.out.marker_genes,
            'MULTI_SAMPLE.final_output',
        )
        
        // Define the parameters for clustering
        def clusteringParams = SC__SCANPY__CLUSTERING_PARAMS( clean(scanpyParams.clustering) )

        // Select a default clustering when in parameter exploration mode
        if(params.tools?.directs && clusteringParams.isParameterExplorationModeOn()) {
            scopeloom = SC__DIRECTS__SELECT_DEFAULT_CLUSTERING( FINALIZE.out.scopeloom )
        } else {
            scopeloom = FINALIZE.out.scopeloom
        }

        // Publishing
        PUBLISH( 
            CLUSTER_IDENTIFICATION.out.marker_genes.map { 
                it -> tuple(it[0], it[1], null)
            },
            params.global.project_name+".multi_sample.final_output",
            "h5ad",
            null,
            clusteringParams.isParameterExplorationModeOn()
        )

        /*******************************************
        * Reporting
        */

        samples = data.map { it -> it[0] }
        UTILS__GENERATE_WORKFLOW_CONFIG_REPORT(
            file(workflow.projectDir + params.utils.workflow_configuration.report_ipynb)
        )

        ipynbs = QC_FILTER.out.report.map {
            it -> tuple(it[0], it[1])
        }.mix(
            samples.combine(UTILS__GENERATE_WORKFLOW_CONFIG_REPORT.out),
            HVG_SELECTION.out.report.map {
                it -> tuple(it[0], it[1])
            },
            DIM_REDUCTION_TSNE_UMAP.out.report.map {
                it -> tuple(it[0], it[1])
            }
        ).groupTuple().map {
            it -> tuple(it[0], *it[1])
        }.join(
            CLUSTER_IDENTIFICATION.out.report,
            by: 0
        )

        if(!clusteringParams.isParameterExplorationModeOn()) {
            ipynbs = ipynbs.map {
                it -> tuple(it[0], it[1..it.size()-2], null)
            }
        } else {
            ipynbs = ipynbs.map {
                it -> tuple(it[0], it[1..it.size()-2], it[it.size()-1])
            }
        }

        SC__SCANPY__MERGE_REPORTS(
            ipynbs,
            "merged_report",
            clusteringParams.isParameterExplorationModeOn()
        )
        SC__SCANPY__REPORT_TO_HTML(SC__SCANPY__MERGE_REPORTS.out)

    emit:
        concatenated_data = scanpyParams?.filter && params.utils?.file_concatenator ? concatenated : Channel.empty()
        final_processed_data = marker_genes
        filteredloom = FINALIZE.out.filteredloom
        scopeloom = scopeloom
        scanpyh5ad = FINALIZE.out.scanpyh5ad

}
