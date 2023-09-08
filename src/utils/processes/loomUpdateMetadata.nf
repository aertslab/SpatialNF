nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")


process SC__UTILS__UPDATE_LOOM_METADATA {
    container params.utils.update_loom_metadata.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__default'

    input:
        tuple val(sampleId), path(f1), path(f2)

    output:
        tuple val(sampleId), path("${sampleId}.SC__UTILS__UPDATE_LOOM_METADATA.loom")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.update_loom_metadata)
        processParams = sampleParams.local

        """
        ${binDir}/update_loom_metadata.py \
            $f1 \
            $f2 \
            "${sampleId}.SC__UTILS__UPDATE_LOOM_METADATA.loom" \
            ${processParams?.axis ? "--axis "+ processParams.axis : "" } \
            ${processParams?.metricKeys ? "--metric-keys "+ processParams.metricKeys : "" } \
            ${processParams?.annotationKeys ? "--annotation-keys "+ processParams.annotationKeys : "" } \
            ${processParams?.embeddingKeys ? "--embedding-keys "+ processParams.embeddingKeys : "" } \
            ${processParams?.clusteringKeys ? "--clustering-keys "+ processParams.clusteringKeys : "" } 
        """

}
