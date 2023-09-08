nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage/bin/" : ""

process KNN_IMPUTATION {

    container params.tools.spage.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1), path(f2)
    output:
        tuple val(sampleId), path("${sampleId}.KNN_IMPUTATION.pkl")
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage)
        processParams = sampleParams.local
        //varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}/knn_imputation.py \
            ${f1} \
            ${f2} \
            ${processParams?.knn ? "--knn "+ processParams.knn : "" } \
            --n-jobs ${task.cpus} \
            ${sampleId}.KNN_IMPUTATION.pkl
        """
}

