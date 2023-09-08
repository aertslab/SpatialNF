
import static groovy.json.JsonOutput.*

nextflow.enable.dsl=2

include { 
    INIT;
} from './src/utils/workflows/utils' params(params)

INIT(params)

include {
    SC__FILE_CONVERTER;
} from './src/utils/processes/utils' params(params)

include {
    getDataChannel;
} from './src/channels/channels' params(params)
include {
    getReferenceDataChannel;
} from './src/channels/channels' params(params)



// run single_sample, output a scope loom file
workflow single_sample {

    include {
        single_sample as SINGLE_SAMPLE;
    } from './workflows/single_sample' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | SINGLE_SAMPLE

    if(params.utils?.publish) {
        PUBLISH_SINGLE_SAMPLE_SCOPE(
            SINGLE_SAMPLE.out.scopeloom,
            "SINGLE_SAMPLE",
            "loom",
            null,
            false
        )
        PUBLISH_SINGLE_SAMPLE_SCANPY(
            SINGLE_SAMPLE.out.scanpyh5ad,
            "SINGLE_SAMPLE",
            "h5ad",
            null,
            false
        )
    }  

}

// run single_sample QC
workflow single_sample_qc {

    include {
        single_sample_qc as SINGLE_SAMPLE_QC;
    } from './src/scanpy/main' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | SINGLE_SAMPLE_QC

    if(params.utils?.publish) {
        PUBLISH(
            SINGLE_SAMPLE_QC.out.filtered,
            "SINGLE_SAMPLE_QC",
            "h5ad",
            null,
            false
        )
    }  

}

workflow multi_sample_qc {

    include {
        multi_sample_qc as MULTI_SAMPLE_QC;
    } from './src/scanpy/main' params(params)
    include {
        PUBLISH;
    } from "./src/utils/workflows/utils" params(params)

    getDataChannel | MULTI_SAMPLE_QC

    if(params.utils?.publish) {
        PUBLISH(
            MULTI_SAMPLE_QC.out.filtered,
            "MULTI_SAMPLE_QC",
            "h5ad",
            null,
            false
        )
    }  

}

workflow multi_sample {

    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)

    getDataChannel | MULTI_SAMPLE
    include {
        PUBLISH as PUBLISH_SCOPE;
        PUBLISH as PUBLISH_SCANPY;
    } from "./src/utils/workflows/utils" params(params)

    if(params.utils?.publish) {
        PUBLISH_SCOPE(
            MULTI_SAMPLE.out.scopeloom,
            "MULTI_SAMPLE",
            "loom",
            null,
            false
        )
        PUBLISH_SCANPY(
            MULTI_SAMPLE.out.scanpyh5ad,
            "MULTI_SAMPLE",
            "h5ad",
            null,
            false
        )
    }

}


// --- SPATIAL WORKFLOWS ---

// run single_sample and spatialDE
workflow single_sample_spatialde {

    include {
        SINGLE_SAMPLE;
    } from "./src/scanpy/workflows/single_sample" params(params)
    include {
        GET_SPATIAL_VARIABLE_GENES as SPATIALDE__GET_SPATIAL_VARIABLE_GENES;
    } from "./src/spatialde/workflows/spatialde" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    
    getDataChannel | SC__FILE_CONVERTER_SPATIAL | SINGLE_SAMPLE 

    input_spatialde = SINGLE_SAMPLE.out.filtered_data.join(SINGLE_SAMPLE.out.final_processed_data)

    SPATIALDE__GET_SPATIAL_VARIABLE_GENES(input_spatialde)

    out = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.out
    scopeout = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.scopeout
    report = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report
    report_aeh = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report_aeh
    
    FILE_CONVERTER_TO_SCOPE(
			scopeout,
			'SINGLE_SAMPLE_SPATIALDE.final_output',
			'mergeToSCopeLoom',			
			SINGLE_SAMPLE.out.filtered_data)
    
    if(params.utils?.publish) {
	PUBLISH_SINGLE_SAMPLE_SCANPY(
            out,
            "SPATIALDE__SPATIAL_VARIABLE_GENES",
            "h5ad",
            null,
            false
        )

	PUBLISH_SINGLE_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "SPATIALDE_scope",
            "loom",
            null,
            false
        )
    }
}

// run multi_sample and spatialDE
workflow multi_sample_spatialde {

    include {
        multi_sample as MULTI_SAMPLE;
    } from "./workflows/multi_sample" params(params)
    include {
        GET_SPATIAL_VARIABLE_GENES as SPATIALDE__GET_SPATIAL_VARIABLE_GENES;
    } from "./src/spatialde/workflows/spatialde" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    
    getDataChannel | MULTI_SAMPLE
    input_spatialde = MULTI_SAMPLE.out.concatenated_data.join(MULTI_SAMPLE.out.final_processed_data)
    
    SPATIALDE__GET_SPATIAL_VARIABLE_GENES(input_spatialde)

    out = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.out
    scopeout = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.scopeout
    report = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report
    report_aeh = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report_aeh
    
    FILE_CONVERTER_TO_SCOPE(
			scopeout,
			'MULTI_SAMPLE_SPATIALDE.final_output',
			'mergeToSCopeLoom',			
			MULTI_SAMPLE.out.concatenated_data)
    
    if(params.utils?.publish) {
	PUBLISH_SINGLE_SAMPLE_SCANPY(
            out,
            "SPATIALDE__SPATIAL_VARIABLE_GENES",
            "h5ad",
            null,
            false
        )

	PUBLISH_SINGLE_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "SPATIALDE_scope",
            "loom",
            null,
            false
        )
    }
}

// run spatialDE starting from filtered raw spatial object h5ad
// currently does not work from loom as clusterings are not properly transferred to object 
workflow spatialde {

    include {
        GET_SPATIAL_VARIABLE_GENES_SIMPLE as SPATIALDE__GET_SPATIAL_VARIABLE_GENES;
    } from "./src/spatialde/workflows/spatialde" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    
    getDataChannel | SC__FILE_CONVERTER_SPATIAL | SPATIALDE__GET_SPATIAL_VARIABLE_GENES

    out = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.out
    scopeout = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.scopeout
    report = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report
    report_aeh = SPATIALDE__GET_SPATIAL_VARIABLE_GENES.out.report_aeh
    
    FILE_CONVERTER_TO_SCOPE(
			scopeout,
			'SPATIALDE.final_output',
			'mergeToSCopeLoomSimple',			
			null)
    
    if(params.utils?.publish) {
	PUBLISH_SINGLE_SAMPLE_SCANPY(
            out,
            "SPATIALDE__SPATIAL_VARIABLE_GENES",
            "h5ad",
            null,
            false
        )

	PUBLISH_SINGLE_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "SPATIALDE_scope",
            "loom",
            null,
            false
        )
    }
}


workflow single_sample_spage {
    include {
        SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        SPAGE__DATA_INTEGRATION;
        SPAGE__LABEL_TRANSFER;
        SPAGE__GENE_IMPUTATION;
    } from "./src/spage/workflows/spage" params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params) 

    data = getDataChannel | SC__FILE_CONVERTER_SPATIAL
    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF

    SINGLE_SAMPLE( data )
    integration_res = SPAGE__DATA_INTEGRATION(
        SINGLE_SAMPLE.out.final_processed_data.map {
            it -> tuple(it[0], it[1])
        },
        ref_data.collectFile()
    )

    if (params.tools?.spage?.label_transfer==true){
        SPAGE__LABEL_TRANSFER( 
            integration_res.data,
            integration_res.knn
        )
        if (params.tools?.spage?.spatial_statistics==true){
            SQUIDPY_ANALYSIS(SPAGE__LABEL_TRANSFER.out)
        } 
        out_h5ad = SPAGE__LABEL_TRANSFER.out
        data = SPAGE__LABEL_TRANSFER.out.combine(ref_data.collectFile())
        
    } else {
        data = integration_res.data
    }

    if (params.tools?.spage?.gene_imputation==true){
        SPAGE__GENE_IMPUTATION(
            data,
            integration_res.knn
        )
        out_h5ad = SPAGE__GENE_IMPUTATION.out
    }

    if (params.tools?.spage?.label_transfer==true || params.tools?.spage?.gene_imputation==true) {
        out_loom = FILE_CONVERTER_TO_SCOPE (
            out_h5ad,
            'SINGLE_SAMPLE__SPAGE.final_output',
            'mergeToSCopeLoom',
            SINGLE_SAMPLE.out.filtered_data
        )

        if(params.utils?.publish) {
            PUBLISH_SINGLE_SAMPLE_SCOPE(
                out_loom,
                "SPAGE__LABEL_TRANSFER",
                "loom",
                null,
                false
            )
            PUBLISH_SINGLE_SAMPLE_SCANPY(
                out_h5ad,
                "SPAGE__LABEL_TRANSFER",
                "h5ad",
                null,
                false
            )
        }
    }    
}


workflow multi_sample_spage {
    include {
        multi_sample as MULTI_SAMPLE;
    } from './workflows/multi_sample' params(params)
    include {
        SPAGE__DATA_INTEGRATION;
        SPAGE__LABEL_TRANSFER;
        SPAGE__GENE_IMPUTATION;
    } from "./src/spage/workflows/spage" params(params)
    include {
        PUBLISH as PUBLISH_MULTI_SAMPLE_SCOPE;
        PUBLISH as PUBLISH_MULTI_SAMPLE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)

    getDataChannel| MULTI_SAMPLE
    out_h5ad = MULTI_SAMPLE.out.scanpyh5ad
    out_loom = MULTI_SAMPLE.out.scopeloom

    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF

    integration_res = SPAGE__DATA_INTEGRATION(
       MULTI_SAMPLE.out.final_processed_data.map {
            it -> tuple(it[0], it[1])
        },
        ref_data.collectFile()
    )
    
    if (params.tools?.spage?.label_transfer==true){ 
        SPAGE__LABEL_TRANSFER(
            integration_res.data,
            integration_res.knn
        )
        if (params.tools?.spage?.spatial_statistics==true){
            SQUIDPY_ANALYSIS(SPAGE__LABEL_TRANSFER.out)
        }
        out_h5ad = SPAGE__LABEL_TRANSFER.out
        data = SPAGE__LABEL_TRANSFER.out.combine(ref_data.collectFile())
    } else {
        data = integration_res.data
    }

    if (params.tools?.spage?.gene_imputation==true){
        SPAGE__GENE_IMPUTATION(
            data,
            integration_res.knn
        )
        out_h5ad = SPAGE__GENE_IMPUTATION.out
    }

    if (params.tools?.spage?.label_transfer==true || params.tools?.spage?.gene_imputation==true) {
        out_loom = FILE_CONVERTER_TO_SCOPE (
            out_h5ad,
            'MULTI_SAMPLE__SPAGE.final_output',
            'mergeToSCopeLoom',
            MULTI_SAMPLE.out.concatenated_data
        )
    } 

    if(params.utils?.publish) {
        PUBLISH_MULTI_SAMPLE_SCOPE(
            out_loom,
            "SPAGE",
            "loom",
            null,
            false
        )
        PUBLISH_MULTI_SAMPLE_SCANPY(
            out_h5ad,
            "SPAGE",
            "h5ad",
            null,
            false
        )
    }
}

workflow spage2vec_spage_label_transfer {
    include {
        SINGLE_SAMPLE;
    } from './src/scanpy/workflows/single_sample' params(params)
    include {
        RUN_SPAGE2VEC;
    } from "./src/spage2vec/workflows/spage2vec" params(params)
    include {
        PSEUDO_CELLS;
    } from './src/spage2vec/processes/pseudo_cells' params(params)
    include {
        SPAGE__DATA_INTEGRATION;
        SPAGE__LABEL_TRANSFER;
        SPAGE__GENE_IMPUTATION;
    } from "./src/spage/workflows/spage" params(params)
    include {
        PUBLISH as PUBLISH_SPAGE2VEC_SPAGE_SCOPE;
        PUBLISH as PUBLISH_SPAGE2VEC_SPAGE_SCANPY;
    } from "./src/utils/workflows/utils" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)

    spage2vec_out = getDataChannel | RUN_SPAGE2VEC
    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF

    pseudocells_out = PSEUDO_CELLS(
        spage2vec_out.map {
            it -> tuple(it[0], it[1])
        }
    )

    SINGLE_SAMPLE( pseudocells_out )

    integration_res = SPAGE__DATA_INTEGRATION(
        SINGLE_SAMPLE.out.final_processed_data.map {
            it -> tuple(it[0], it[1])
        },
        ref_data.collectFile()
    )

    if (params.tools?.spage?.label_transfer==true){
        SPAGE__LABEL_TRANSFER(
            integration_res.data,
            integration_res.knn
        )
        if (params.tools?.spage?.spatial_statistics==true){
            SQUIDPY_ANALYSIS(SPAGE__LABEL_TRANSFER.out)
        }
        out_h5ad = SPAGE__LABEL_TRANSFER.out
        data = SPAGE__LABEL_TRANSFER.out.combine(ref_data.collectFile())

    } else {
        data = integration_res.data
    }

    if (params.tools?.spage?.gene_imputation==true){
        SPAGE__GENE_IMPUTATION(
            data,
            integration_res.knn
        )
        out_h5ad = SPAGE__GENE_IMPUTATION.out
    }

    if (params.tools?.spage?.label_transfer==true || params.tools?.spage?.gene_imputation==true) {
        out_loom = FILE_CONVERTER_TO_SCOPE (
            out_h5ad,
            'SPAGE2VEC__SPAGE.final_output',
            'mergeToSCopeLoom',
            SINGLE_SAMPLE.out.filtered_data           
        )
    }

    if(params.utils?.publish) {
        PUBLISH_SPAGE2VEC_SPAGE_SCOPE(
            out_loom,
            "SPAGE2VEC_SPAGE",
            "loom",
            null,
            false
        )
        PUBLISH_SPAGE2VEC_SPAGE_SCANPY(
            out_h5ad,
            "SPAGE2VEC_SPAGE",
            "h5ad",
            null,
            false
        )
    }
}

workflow single_sample_tangram {
    include {
        SINGLE_SAMPLE;
    } from "./src/scanpy/workflows/single_sample" params(params)
    include {
        PROJECT_CELLTYPES as TANGRAM__MAP_CELLTYPES;
    } from "./src/tangram/workflows/tangram" params(params)
    include {
    	TANGRAM__PREPARE_SCRNA as TANGRAM__PERPARE_REF;
    } from "./src/tangram/processes/run_tangram.nf" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
	PUBLISH as PUBLISH_SINGLE_SAMPLE_MAPPED;
	PUBLISH as PUBLISH_SINGLE_SAMPLE_MAPPING;
    } from "./src/utils/workflows/utils" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)


    data = getDataChannel | SC__FILE_CONVERTER_SPATIAL
    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF | TANGRAM__PERPARE_REF

    SINGLE_SAMPLE (data)

    input_spatial = SINGLE_SAMPLE.out.final_processed_data.combine(SINGLE_SAMPLE.out.filtered_data, by: 0)
    TANGRAM__MAP_CELLTYPES( input_spatial.combine(ref_data) )

    FILE_CONVERTER_TO_SCOPE(
			TANGRAM__MAP_CELLTYPES.out.mapped,
			'SINGLE_SAMPLE_TANGRAM.final_output',
			'mergeToSCopeLoomSimple',
			null)

    if (params.tools?.tangram?.squidpy_statistics==true){
            SQUIDPY_ANALYSIS(TANGRAM__MAP_CELLTYPES.out.mapped)
        }

    if(params.utils?.publish) {
        PUBLISH_SINGLE_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "TANGRAM_CELLTYPES_scope",
            "loom",
            null,
            false
        )
	PUBLISH_SINGLE_SAMPLE_MAPPING(
	    TANGRAM__MAP_CELLTYPES.out.mapping,
            "TANGRAM_MAPPING",
            "h5ad",
            null,
            false
        )
	PUBLISH_SINGLE_SAMPLE_MAPPED(
	    TANGRAM__MAP_CELLTYPES.out.mapped,
            "TANGRAM_CELLTYPES",
            "h5ad",
            null,
            false
        )
    }
}


workflow tangram {
    include {
        PROJECT_CELLTYPES_SIMPLE as TANGRAM__MAP_CELLTYPES;
    } from "./src/tangram/workflows/tangram" params(params)
    include {
    	TANGRAM__PREPARE_SCRNA as TANGRAM__PERPARE_REF;
    } from "./src/tangram/processes/run_tangram.nf" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_SINGLE_SAMPLE_SCOPE;
	PUBLISH as PUBLISH_SINGLE_SAMPLE_MAPPED;
	PUBLISH as PUBLISH_SINGLE_SAMPLE_MAPPING;
    } from "./src/utils/workflows/utils" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_SPATIAL;
    } from './src/utils/processes/utils' params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)


    data = getDataChannel | SC__FILE_CONVERTER_SPATIAL
    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF | TANGRAM__PERPARE_REF

    TANGRAM__MAP_CELLTYPES( data.combine(ref_data) )

    FILE_CONVERTER_TO_SCOPE(
			TANGRAM__MAP_CELLTYPES.out.mapped,
			'SINGLE_SAMPLE_TANGRAM.final_output',
			'mergeToSCopeLoomSimple',
			null)

    if (params.tools?.tangram?.squidpy_statistics==true){
            SQUIDPY_ANALYSIS(TANGRAM__MAP_CELLTYPES.out.mapped)
        }

    if(params.utils?.publish) {
        PUBLISH_SINGLE_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "TANGRAM_CELLTYPES_scope",
            "loom",
            null,
            false
        )
	PUBLISH_SINGLE_SAMPLE_MAPPING(
	    TANGRAM__MAP_CELLTYPES.out.mapping,
            "TANGRAM_MAPPING",
            "h5ad",
            null,
            false
        )
	PUBLISH_SINGLE_SAMPLE_MAPPED(
	    TANGRAM__MAP_CELLTYPES.out.mapped,
            "TANGRAM_CELLTYPES",
            "h5ad",
            null,
            false
        )
    }
}

workflow multi_sample_tangram {
    include {
        multi_sample as MULTI_SAMPLE;
    } from "./workflows/multi_sample" params(params)
    include {
        PROJECT_CELLTYPES as TANGRAM__MAP_CELLTYPES;
    } from "./src/tangram/workflows/tangram" params(params)
    include {
    	TANGRAM__PREPARE_SCRNA as TANGRAM__PERPARE_REF;
    } from "./src/tangram/processes/run_tangram.nf" params(params)
    include {
        FILE_CONVERTER as FILE_CONVERTER_TO_SCOPE;
    } from "./src/utils/workflows/fileConverter" params(params)
    include {
        PUBLISH as PUBLISH_MULTI_SAMPLE_SCOPE;
	PUBLISH as PUBLISH_MULTI_SAMPLE_MAPPED;
	PUBLISH as PUBLISH_MULTI_SAMPLE_MAPPING;
    } from "./src/utils/workflows/utils" params(params)
    include {
    	SC__FILE_CONVERTER as SC__FILE_CONVERTER_REF;
    } from './src/utils/processes/utils' params(params)
    include {
        SQUIDPY_ANALYSIS;
    } from './src/squidpy/workflows/squidpy_analysis' params(params)

    getDataChannel | MULTI_SAMPLE
    ref_data = getReferenceDataChannel | SC__FILE_CONVERTER_REF | TANGRAM__PERPARE_REF

    input_spatial = MULTI_SAMPLE.out.final_processed_data.combine(MULTI_SAMPLE.out.concatenated_data, by: 0)
    
    TANGRAM__MAP_CELLTYPES( input_spatial.combine(ref_data) )

    FILE_CONVERTER_TO_SCOPE(
			TANGRAM__MAP_CELLTYPES.out.mapped,
			'MULTI_SAMPLE_TANGRAM.final_output',
			'mergeToSCopeLoomSimple',
			null)

    if (params.tools?.tangram?.squidpy_statistics==true){
            SQUIDPY_ANALYSIS(TANGRAM__MAP_CELLTYPES.out.mapped)
        }

    if(params.utils?.publish) {
        PUBLISH_MULTI_SAMPLE_SCOPE(
	    FILE_CONVERTER_TO_SCOPE.out,
            "TANGRAM_CELLTYPES_scope",
            "loom",
            null,
            false
        )
	PUBLISH_MULTI_SAMPLE_MAPPING(
	    TANGRAM__MAP_CELLTYPES.out.mapping,
            "TANGRAM_MAPPING",
            "h5ad",
            null,
            false
        )
	PUBLISH_MULTI_SAMPLE_MAPPED(
	    TANGRAM__MAP_CELLTYPES.out.mapped,
            "TANGRAM_CELLTYPES",
            "h5ad",
            null,
            false
        )
    }
}
