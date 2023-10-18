nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage2vec/bin/" : ""

process PSEUDO_CELLS {

    container params.tools.spage2vec.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.PSEUDO_CELLS.h5ad")
    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage2vec)
        processParams = sampleParams.local
        """
        ${binDir}/pseudo_cells.py \
            ${f} \
            ${processParams?.n_neigh ? "--n_neigh "+ processParams.n_neigh : "" } \
            ${processParams?.n_cpus ? "--n_cpus "+ processParams.n_cpus : "" } \
            ${sampleId}.PSEUDO_CELLS.h5ad \
        """
}
