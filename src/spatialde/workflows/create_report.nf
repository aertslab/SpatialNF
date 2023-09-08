nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SPATIALDE__SCANPY__GENERATE_REPORT as SPATIALDE__GENERATE_REPORT;
    SPATIALDE__SCANPY__REPORT_TO_HTML as SPATIALDE__REPORT_TO_HTML; 
} from '../processes/reports.nf' params(params)



workflow GENERATE_REPORT {

    take:
        pipelineStep
        data // anndata
        ipynb

    main:
        def reportTitle = "spatialDE_" + pipelineStep.toLowerCase() + "_report"

	report_notebook = SPATIALDE__GENERATE_REPORT(
                        ipynb,
                        // expects (sample_id, adata)
                        data,
                        reportTitle
                    ).map {
                        it -> tuple(it[0], it[1], null)
                    }
		    
        SPATIALDE__REPORT_TO_HTML(report_notebook.map {
            it -> tuple(it[0], it[1])
        })

    emit:
        report_notebook
}
