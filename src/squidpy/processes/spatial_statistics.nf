nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/squidpy/bin/" : ""

process SPATIAL_STATISTICS {

    container params.tools.squidpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1)
    output:
        tuple val(sampleId), path("${sampleId}.SPATIAL_STATISTICS.${processParams.spatial_statistics.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.squidpy)
        processParams = sampleParams.local
        """
        ${binDir}/spatial_statistics.py \
            ${f1} \
            ${processParams?.spatial_statistics.cluster_key ? "--cluster_key "+ processParams.spatial_statistics.cluster_key : "" } \
            ${processParams?.spatial_statistics.connectivity_key ? "--connectivity_key "+ processParams.spatial_statistics.connectivity_key : "" } \
            ${processParams?.spatial_statistics.n_perms ? "--n_perms "+ processParams.spatial_statistics.n_perms : "" } \
            ${processParams?.spatial_statistics.ripley_mod ? "--ripley_mod "+ processParams.spatial_statistics.ripley_mod : "" } \
            ${processParams?.spatial_statistics.n_steps ? "--n_steps "+ processParams.spatial_statistics.n_steps : "" } \
            ${processParams?.spatial_statistics.n_neigh ? "--n_neigh "+ processParams.spatial_statistics.n_neigh : "" } \
            ${processParams?.spatial_statistics.n_simulations ? "--n_simulations "+ processParams.spatial_statistics.n_simulations : "" } \
            ${processParams?.spatial_statistics.n_observations ? "--n_observations "+ processParams.spatial_statistics.n_observations : "" } \
            ${processParams?.spatial_statistics.max_dist ? "--max_dist "+ processParams.spatial_statistics.max_dist : "" } \
            ${params.global?.seed ? "--seed "+ params.global.seed : ""} \
            --n_jobs ${task.cpus} \
            ${sampleId}.SPATIAL_STATISTICS.${processParams.spatial_statistics.off} \
        """
}

