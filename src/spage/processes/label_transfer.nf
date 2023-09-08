nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage/bin/" : ""

process LABEL_TRANSFER {

    container params.tools.spage.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1), path(f2), path(f3)
    output:
        tuple val(sampleId), path("${sampleId}.LABEL_TRANSFER.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage)
        processParams = sampleParams.local
        //varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}/label_transfer.py \
            ${f1} \
            ${f2} \
            ${f3} \
            ${processParams?.label_key ? "--label-key "+ processParams.label_key : "" } \
            ${sampleId}.LABEL_TRANSFER.${processParams.off} \
        """
}

