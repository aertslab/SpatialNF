nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_TANGRAM_CELLTYPES;
} from "../../utils/workflows/utils.nf" params(params)
include {
    PUBLISH as PUBLISH_H5AD_TANGRAM_MAPPING;
} from "../../utils/workflows/utils.nf" params(params)

include {
    TANGRAM__PREPARE_SCRNA;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__PREPARE_SPATIAL;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__PREPARE_SPATIAL_SIMPLE;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__COMPUTE_MAPPING;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__PROJECT_CELLTYPES;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__PROJECT_EXPRESSION;
} from '../processes/run_tangram.nf' params(params)
include {
    TANGRAM__GENERATE_DUAL_INPUT_REPORT; 
} from '../processes/reports.nf' params(params)
include {
    TANGRAM__REPORT_TO_HTML; 
} from '../processes/reports.nf' params(params)


//////////////////////////////////////////////////////


workflow PROJECT_CELLTYPES {

    take:
	input

    main:

	input.multiMap { it ->
		       data: tuple(it[0], it[1])
		       filtered: tuple(it[0], it[3])
		       ref_data: tuple(it[4], it[5])
		       }
		       .set{ combData }

        data = combData.data
	filtered = combData.filtered
	ref = combData.ref_data
	
	spatial = TANGRAM__PREPARE_SPATIAL( data, filtered )
	mapping = TANGRAM__COMPUTE_MAPPING(spatial, ref)
        mapped  = TANGRAM__PROJECT_CELLTYPES( spatial, ref, mapping)

	if(params.tools?.tangram?.project_gex==true){
		gex = TANGRAM__PROJECT_EXPRESSION( mapped, ref, mapping)
		mapped = gex
	}
	
        PUBLISH_H5AD_TANGRAM_CELLTYPES(
		mapped.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "TANGRAM.celltypes_output",
            "h5ad",
            "tangram",
            false
        )
	PUBLISH_H5AD_TANGRAM_MAPPING(
		mapping.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "TANGRAM.mapping_output",
            "h5ad",
            "tangram",
            false
        )

	report = Channel.empty()
	if(params.tools?.tangram?.report_ipynb)
	{
	    report = TANGRAM__GENERATE_DUAL_INPUT_REPORT(
		file(workflow.projectDir + params.tools.tangram.report_ipynb),
                mapping.join(mapped).map { 
                    it -> tuple(*it[0..(it.size()-1)], null)
                },                
                'TANGRAM_report',
                false
            )
	    TANGRAM__REPORT_TO_HTML(report.map {
            	it -> tuple(it[0], it[1])
            })
	}
	
    emit:
	mapped
	mapping
	report
}


workflow PROJECT_CELLTYPES_SIMPLE {

    take:
	input

    main:

	input.multiMap { it ->
		       data: tuple(it[0], it[1])
		       ref_data: tuple(it[2], it[3])
		       }
		       .set{ combData }

        data = combData.data
	ref = combData.ref_data
	
	spatial = TANGRAM__PREPARE_SPATIAL_SIMPLE( data )
	mapping = TANGRAM__COMPUTE_MAPPING(spatial, ref)
        mapped  = TANGRAM__PROJECT_CELLTYPES( spatial, ref, mapping)

	if(params.tools?.tangram?.project_gex==true){
		gex = TANGRAM__PROJECT_EXPRESSION( mapped, ref, mapping)
		mapped = gex
	}
	
        PUBLISH_H5AD_TANGRAM_CELLTYPES(
		mapped.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "TANGRAM.celltypes_output",
            "h5ad",
            "tangram",
            false
        )
	PUBLISH_H5AD_TANGRAM_MAPPING(
		mapping.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "TANGRAM.mapping_output",
            "h5ad",
            "tangram",
            false
        )

	report = Channel.empty()
	if(params.tools?.tangram?.report_ipynb)
	{
	    report = TANGRAM__GENERATE_DUAL_INPUT_REPORT(
		file(workflow.projectDir + params.tools.tangram.report_ipynb),
                mapping.join(mapped).map { 
                    it -> tuple(*it[0..(it.size()-1)], null)
                },                
                'TANGRAM_report',
                false
            )
	    TANGRAM__REPORT_TO_HTML(report.map {
            	it -> tuple(it[0], it[1])
            })
	}
	
    emit:
	mapped
	mapping
	report
}
