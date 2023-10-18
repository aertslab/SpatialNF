nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_SPAGE2VEC;
} from "../../utils/workflows/utils.nf" params(params)
include {
    READ_COORD_CONVERTER;
} from '../processes/read_coord_converter.nf' params(params)
include {
    READ_COORD_CONCATENATOR;
} from '../processes/read_coord_concatenator.nf' params(params)
include {
    SPAGE2VEC__TRAINING;
} from '../processes/spage2vec_train.nf' params(params)
include {
    getBaseName;
} from '../../utils/processes/files.nf'



workflow RUN__MULTISAMPLE_SPAGE2VEC {
    take:
        data
    main:
        out = READ_COORD_CONVERTER(data)
        concatenated = READ_COORD_CONCATENATOR(
                out.map {
                    it -> it[1]
                }.toSortedList(
                    { a, b -> getBaseName(a, "READ_COORD") <=> getBaseName(b, "READ_COORD") }
                )
        )

        spage2vec_out = SPAGE2VEC__TRAINING(concatenated)
    emit:
        spage2vec_out
}

workflow RUN__SINGLESAMPLE_SPAGE2VEC {
    take:
        data
    main:
        out = READ_COORD_CONVERTER(data)
        spage2vec_out = SPAGE2VEC__TRAINING(out)
    emit:
        spage2vec_out
}
