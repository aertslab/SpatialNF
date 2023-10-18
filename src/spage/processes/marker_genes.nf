nextflow.enable.dsl=2

import java.nio.file.Paths

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin" : Paths.get(workflow.scriptFile.getParent().toString(), "bin")

process MARKER_GENES {
        container params.tools.scanpy.container
        publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
        label 'compute_resources__mem'

        input:
                // Expects 
                // - normalizedTransformedData to be an AnnData containing log normalised data
        tuple \
                        val(sampleId), \
                        path(clusteredData)

        output:
        tuple val(sampleId), path("${sampleId}.SPAGE__MARKER_GENES.${processParams.off}")

        script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spage.marker_genes)
                processParams = sampleParams.local
                """
                ${binDir}/cluster/sc_marker_genes.py \
                        ${(processParams.containsKey('method')) ? '--method ' + processParams.method : ''} \
                        ${(processParams.containsKey('groupby')) ? '--groupby ' + processParams.groupby : ''} \
                        ${(processParams.containsKey('ngenes')) ? '--ngenes ' + processParams.ngenes : ''} \
                        $clusteredData \
                        $clusteredData \
                        "${sampleId}.SPAGE__MARKER_GENES.${processParams.off}"
                """

}
