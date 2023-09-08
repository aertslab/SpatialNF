nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

include {
    SPATIAL__SQUIDPY__GENERATE_REPORT;
    SPATIAL__SQUIDPY__REPORT_TO_HTML;
} from '../../squidpy/processes/reports.nf' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow GENERATE_REPORT {

    take:
        data // anndata
        ipynb

    main:
        def reportTitle = "Squidpy_spatial_statistics_report"
        report_notebook = SPATIAL__SQUIDPY__GENERATE_REPORT(
            ipynb,
            // expects (sample_id, adata)
            data,
            reportTitle
            ).map {
                it -> tuple(it[0], it[1], null)
            }

        SPATIAL__SQUIDPY__REPORT_TO_HTML(report_notebook.map {
            it -> tuple(it[0], it[1])
        })
    emit:
        report_notebook

}

