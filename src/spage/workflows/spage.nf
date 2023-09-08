nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_LABEL_TRANSFER;
    PUBLISH as PUBLISH_H5AD_GENE_IMPUTATION;
} from "../../utils/workflows/utils.nf" params(params)

include {
    DATA_INTEGRATION;
} from '../processes/data_integration.nf' params(params)

include {
    KNN_IMPUTATION;
} from '../processes/knn_imputation.nf' params(params)

include {
    LABEL_TRANSFER;
} from '../processes/label_transfer.nf' params(params)

include {
    GENE_IMPUTATION;
} from '../processes/gene_imputation.nf' params(params)

//////////////////////////////////////////////////////

workflow SPAGE__DATA_INTEGRATION {

    take:
        spatial_data
        sc_data
    main:
        data = DATA_INTEGRATION( spatial_data, sc_data )
        knn = KNN_IMPUTATION( data )
    emit:
        data
        knn
}

workflow SPAGE__LABEL_TRANSFER {

    take:
        data
        knn_data
    main:
        LABEL_TRANSFER( data.join(knn_data) )

        PUBLISH_H5AD_LABEL_TRANSFER(
            LABEL_TRANSFER.out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPAGE.label_transfer_output",
            "h5ad",
            "spage",
            false
        )
    emit:
        LABEL_TRANSFER.out
}

workflow SPAGE__GENE_IMPUTATION {

    take:
        data
        knn_data
    main:
        GENE_IMPUTATION( data.join(knn_data) )

        PUBLISH_H5AD_GENE_IMPUTATION(
            GENE_IMPUTATION.out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPAGE.gene_imputation_output",
            "h5ad",
            "spage",
            false
        )

    emit:
        GENE_IMPUTATION.out
}
