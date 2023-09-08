nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_SPATIAL_VARIABLE_GENES;
} from "../../utils/workflows/utils.nf" params(params)


include {
    SPATIALDE__VARIABLE_GENES;
} from '../processes/run_spatialde.nf' params(params)
include {
    SPATIALDE__SPATIAL_PATTERNS;
} from '../processes/run_spatialde.nf' params(params)
include {
        SPATIALDE__ADD_SPATIAL_PATTERNS;
    } from "../processes/run_spatialde" params(params)
include {
    GENERATE_REPORT as GENERATE_REPORT_VARIABLE_GENES
    GENERATE_REPORT as GENERATE_REPORT_AEH;	
} from "./create_report.nf" params(params)

//////////////////////////////////////////////////////

workflow GET_SPATIAL_VARIABLE_GENES {

    take:
        input
	// input: single sample workflow: 'filtered_data' joined with 'final_processed_data'
	
    main:
	def spatialdeParams = params.tools.spatialde

	input.multiMap { it ->
                       filtered_data: tuple(it[0], it[1])
                       final_processed_data: tuple(it[0], it[2])
                       }
                       .set{ data }	

	out = SPATIALDE__VARIABLE_GENES( data.filtered_data )
	
        PUBLISH_H5AD_SPATIAL_VARIABLE_GENES(
            out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPATIALDE.variable_genes_output",
            "h5ad",
            "spatialDE",
            false
        )

	// init empty reports to ensure return values
	report = null
	report_aeh = null
	
	if(spatialdeParams?.report_ipynb)
        {
	    report = GENERATE_REPORT_VARIABLE_GENES("variable_genes", out, file(workflow.projectDir + spatialdeParams.report_ipynb) )
        }

        if(spatialdeParams?.run_aeh)
        {
	    out = SPATIALDE__SPATIAL_PATTERNS(out)
	    
	    if(spatialdeParams?.report_aeh_ipynb)
    	    {
		report_aeh = GENERATE_REPORT_AEH("aeh", out, file(workflow.projectDir + spatialdeParams.report_aeh_ipynb))
    	    }
	
	    scopeout = SPATIALDE__ADD_SPATIAL_PATTERNS(out, data.final_processed_data  )	
        } else {
             scopeout = data.final_processed_data 
        }

    emit:
        out
	scopeout
	report
	report_aeh
}


workflow GET_SPATIAL_VARIABLE_GENES_SIMPLE {

    take:
        input
	// input: spatial data object
	
    main:
	def spatialdeParams = params.tools.spatialde

	input.multiMap { it ->
                       filtered_data: tuple(it[0], it[1])
                       }
                       .set{ data }	

	out = SPATIALDE__VARIABLE_GENES( data.filtered_data )
	
        PUBLISH_H5AD_SPATIAL_VARIABLE_GENES(
            out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPATIALDE.variable_genes_output",
            "h5ad",
            "spatialDE",
            false
        )

	// init empty reports to ensure return values
	report = null
	report_aeh = null
	
	if(spatialdeParams?.report_ipynb)
        {
	    report = GENERATE_REPORT_VARIABLE_GENES("variable_genes", out, file(workflow.projectDir + spatialdeParams.report_ipynb) )
        }

        if(spatialdeParams?.run_aeh)
        {
	    out = SPATIALDE__SPATIAL_PATTERNS(out)
	    
	    if(spatialdeParams?.report_aeh_ipynb)
    	    {
		report_aeh = GENERATE_REPORT_AEH("aeh", out, file(workflow.projectDir + spatialdeParams.report_aeh_ipynb))
    	    }
	
	    scopeout = SPATIALDE__ADD_SPATIAL_PATTERNS(out, data.filtered_data  )	
        } else {
             scopeout = data.filtered_data
        }

    emit:
        out
	scopeout
	report
	report_aeh
}
