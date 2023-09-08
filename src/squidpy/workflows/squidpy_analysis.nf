nextflow.enable.dsl=2

//////////////////////////////////////////////////////
include {
    SPATIAL_GRAPH;
} from '../processes/neighborhood_graph.nf' params(params)
include {
    SPATIAL_STATISTICS;
} from '../processes/spatial_statistics.nf' params(params)
include {
    LR_ANALYSIS;
} from '../processes/lr_analysis.nf' params(params)
// Reporting
include {
    GENERATE_REPORT;
} from './create_report.nf' params(params + params.global)
include {
    PUBLISH as PUBLISH_FINAL_SQUIDPY_OUTPUT;
} from "../../utils/workflows/utils.nf" params(params)

//////////////////////////////////////////////////////

workflow SQUIDPY_ANALYSIS {

    take:
        data
    main:
        // To avoid Variable `params` already defined in the process scope
        def squidpyParams = params.tools.squidpy

        out = SPATIAL_GRAPH( data )
        out_statistics = SPATIAL_STATISTICS( out )
        out = LR_ANALYSIS( out_statistics )

	// Publish final h5ad
        //PUBLISH_FINAL_SQUIDPY_OUTPUT(
        //    out.map { it -> tuple(it[0], it[1]) },
        //    'SQUIDPY_ANALYSIS.output',
        //    "h5ad",
        //     null,
        //    false
        //)
        squidpy_report = GENERATE_REPORT(
            out_statistics.map { it -> tuple(it[0], it[1]) },
            file(workflow.projectDir + params.tools.squidpy.spatial_statistics.report_ipynb),
        )
        
    emit:
        out
        squidpy_report
}
