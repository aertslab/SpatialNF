nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage/bin/" : ""

process DATA_INTEGRATION {

    container params.tools.spage.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1)
        path(f2)
    output:
        tuple val(sampleId), \
            path("${sampleId}.DATA_INTEGRATION_SPATIAL.${processParams.off}"), \
            path("${sampleId}.DATA_INTEGRATION_SC.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage)
        processParams = sampleParams.local
        //varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}/integrate.py \
            ${f1} \
            ${f2} \
            ${processParams?.npvs ? "--npvs "+ processParams.npvs : "" } \
            ${processParams?.th ? "--th "+ processParams.th : "" } \
            ${processParams?.normalize_sc ? "--norm "+ processParams.normalize_sc : "" } \
            ${sampleId}.DATA_INTEGRATION_SPATIAL.${processParams.off} \
            ${sampleId}.DATA_INTEGRATION_SC.${processParams.off}
        """
}

