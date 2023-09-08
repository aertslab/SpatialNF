nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:

// Imports
include {
    SC__UTILS__UPDATE_LOOM_METADATA;
} from './../processes/loomUpdateMetadata' params(params)

//////////////////////////////////////////////////////
//  Define the workflow 

workflow UPDATE_LOOM_METADATA {

    take:
        // Expects (sampleId, data)
        data

    main:
        out = SC__UTILS__UPDATE_LOOM_METADATA( data )
    emit:
        out

}
