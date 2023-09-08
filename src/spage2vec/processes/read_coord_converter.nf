nextflow.enable.dsl=2
binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spage2vec/bin/" : ""

process READ_COORD_CONVERTER {

    cache 'deep'
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    container params.tools.spage2vec.container
    label 'compute_resources__mem'

    input:
        tuple val(sampleId), path(f), val(NULL), val(NULL), val(NULL)

    output:
        tuple val(sampleId), path("${sampleId}.READ_COORD__OUTPUT.csv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.file_converter)
        processParams = sampleParams.local

        """
        ${binDir}/read_coord_converter.py \
            ${f} \
            --sample_id ${sampleId} \
            ${processParams?.x_col_name ? "--x_col_name "+ processParams.x_col_name : "" } \
            ${processParams?.y_col_name ? "--y_col_name "+ processParams.y_col_name : "" } \
            ${processParams?.z_col_name ? "--z_col_name "+ processParams.z_col_name : "" } \
            ${processParams?.gene_col_name ? "--gene_col_name "+ processParams.gene_col_name : "" } \
            ${processParams?.delimiter ? "--delimiter "+ processParams.delimiter : "" } \
            ${sampleId}.READ_COORD__OUTPUT.csv
        """
}
