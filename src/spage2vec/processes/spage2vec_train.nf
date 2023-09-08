nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage2vec/bin/" : ""

process SPAGE2VEC__TRAINING {

    container params.tools.spage2vec.container
    publishDir "${params.global.outdir}/spage2vec/model", mode: 'symlink'
    label 'compute_resources__gpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), \
            path("${sampleId}.SPAGE2VEC__OUTPUT.h5ad"), \
            file("${sampleId}.SPAGE2VEC__MODEL.h5")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage2vec)
        processParams = sampleParams.local
        """
        ${binDir}/spage2vec_train.py \
            ${f} \
            ${processParams?.feature_key ? "--feature_key "+ processParams.feature_key : "" } \
            ${processParams?.spatial_key ? "--spatial_key "+ processParams.spatial_key : "" } \
            ${processParams?.radius ? "--radius "+ processParams.radius : "" } \
            ${processParams?.neigh_perc ? "--neigh_perc "+ processParams.neigh_perc : "" } \
            ${processParams?.delaunay ? "--delaunay "+ processParams.delaunay : "" } \
            ${processParams?.min_conn_comp_size ? "--min_conn_comp_size "+ processParams.min_conn_comp_size : "" } \
            ${processParams?.num_of_walks ? "--num_of_walks "+ processParams.num_of_walks : "" } \
            ${processParams?.walk_length ? "--walk_length "+ processParams.walk_length : "" } \
            ${processParams?.batch_size ? "--batch_size "+ processParams.batch_size : "" } \
            ${processParams?.num_epochs ? "--num_epochs "+ processParams.num_epochs : "" } \
            ${processParams?.k1 ? "--k1 "+ processParams.k1 : "" } \
            ${processParams?.k2 ? "--k2 "+ processParams.k2 : "" } \
            ${processParams?.l1 ? "--l1 "+ processParams.l1 : "" } \
            ${processParams?.l2 ? "--l2 "+ processParams.l2 : "" } \
            ${processParams?.device ? "--gpu_id "+ processParams.device : "" } \
            --n_cpus ${task.cpus} \
            --seed ${params.global.seed} \
            ${sampleId}.SPAGE2VEC__OUTPUT.h5ad \
            ${sampleId}.SPAGE2VEC__MODEL.h5
        """
}

