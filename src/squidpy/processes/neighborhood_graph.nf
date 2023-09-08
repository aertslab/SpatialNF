nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/squidpy/bin/" : ""

process SPATIAL_GRAPH {

    container params.tools.squidpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f1)
    output:
        tuple val(sampleId), path("${sampleId}.SPATIAL_GRAPH.${processParams.spatial_graph.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.squidpy)
        processParams = sampleParams.local
        """
        ${binDir}/neighborhood_graph.py \
            ${f1} \
            ${processParams?.spatial_graph.spatial_key ? "--spatial_key "+ processParams.spatial_graph.spatial_key : "" } \
            ${processParams?.spatial_graph.coord_type ? "--coord_type "+ processParams.spatial_graph.coord_type : "" } \
            ${processParams?.spatial_graph.n_neighs ? "--n_neighs "+ processParams.spatial_graph.n_neighs : "" } \
            ${processParams?.spatial_graph.radius ? "--radius "+ processParams.spatial_graph.radius : "" } \
            ${processParams?.spatial_graph.delauny ? "--delauny "+ processParams.spatial_graph.delauny : "" } \
            ${processParams?.spatial_graph.n_rings ? "--n_rings "+ processParams.spatial_graph.n_rings : "" } \
            ${processParams?.spatial_graph.set_diag ? "--set_diag "+ processParams.spatial_graph.set_diag : "" } \
            ${sampleId}.SPATIAL_GRAPH.${processParams.spatial_graph.off} \
        """
}


