nextflow.enable.dsl=2
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage2vec/bin/" : ""

process READ_COORD_CONCATENATOR {

    cache 'deep'
    container params.tools.spage2vec.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        file("*")

    output:
        tuple val(params.global.project_name), path("${params.global.project_name}.READ_COORD_CONCATENATOR.csv")

    script:
        """
        ${binDir}/read_coord_concatenator.py * "${params.global.project_name}.READ_COORD_CONCATENATOR.csv" 
        """

}
