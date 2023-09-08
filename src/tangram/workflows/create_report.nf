nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    TANGRAM__GENERATE_DUAL_INPUT_REPORT; 
} from '../processes/reports.nf' params(params)
include {
    SC__SCANPY__REPORT_TO_HTML as TANGRAM__REPORT_TO_HTML; 
} from '../../scanpy/processes/reports.nf' params(params)



workflow GENERATE_REPORT {

    take:
        pipelineStep
        data1 // anndata1
	data2 // anndata2
        ipynb

    main:
        def reportTitle = "tangram_" + pipelineStep.toLowerCase() + "_report"

        report_notebook = TANGRAM__GENERATE_DUAL_INPUT_REPORT(
                        ipynb,
                        data1,
			data2,
			
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
