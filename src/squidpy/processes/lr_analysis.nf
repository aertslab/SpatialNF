nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/squidpy/bin/" : ""

process LR_ANALYSIS {

    container params.tools.squidpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1)
    output:
        tuple val(sampleId), path("${sampleId}.LR_ANALYSIS.${processParams.lr_analysis.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.squidpy)
        processParams = sampleParams.local
        """
        ${binDir}/lr_analysis.py \
            ${f1} \
            ${processParams?.lr_analysis.cluster_key ? "--cluster_key "+ processParams.lr_analysis.cluster_key : "" } \
            ${processParams?.lr_analysis.threshold ? "--threshold "+ processParams.lr_analysis.threshold : "" } \
            ${processParams?.lr_analysis.n_perms ? "--n_perms "+ processParams.lr_analysis.n_perms : "" } \
            ${processParams?.lr_analysis.complex_policy ? "--complex_policy "+ processParams.lr_analysis.complex_policy : "" } \
            ${processParams?.lr_analysis.corr_method ? "--corr_method "+ processParams.lr_analysis.corr_method : "" } \
            ${processParams?.lr_analysis.corr_axis ? "--corr_axis "+ processParams.lr_analysis.corr_axis : "" } \
            ${processParams?.lr_analysis.corr_alpha ? "--corr_alpha "+ processParams.lr_analysis.corr_alpha : "" } \
	    ${processParams?.lr_analysis.not_use_raw ? "--not_use_raw "+ processParams.lr_analysis.not_use_raw : "" } \
            ${params.global?.seed ? "--seed "+ params.global.seed : "" } \
            --n_jobs ${task.cpus} \
            ${sampleId}.LR_ANALYSIS.${processParams.lr_analysis.off} \
        """
}
