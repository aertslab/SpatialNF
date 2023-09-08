nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spatialde/bin/" : ""

process SPATIALDE__VARIABLE_GENES {
    
    container params.tools.spatialde.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SPATIALDE__VARIABLE_GENES.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spatialde)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/run_spatialde.py \
		$f \
      		"${sampleId}.SPATIALDE__VARIABLE_GENES.${processParams.off}" \
		${processParams?.min_cells ? "--min-cells " + processParams.min_cells  : ""} \
		${processParams?.thr_qval ? "--thr-qval " + processParams.thr_qval  : ""} \
		${processParams?.method_pval_correction ? "--method-pval-correction " + processParams.method_pval_correction  : ""} \
		${processParams?.min_fsv ? "--min-fsv " + processParams.min_fsv  : ""} \
		${processParams?.normalize_naivede ? "--normalize-naivede " + processParams.normalize_naivede  : ""} \
		--ncpu ${task.cpus}
      	"""
}


process SPATIALDE__SPATIAL_PATTERNS {
    container params.tools.spatialde.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SPATIALDE__SPATIAL_PATTERNS.${processParams.off}")
	
    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spatialde)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/run_aeh.py \
		$f \
      		${sampleId}.SPATIALDE__SPATIAL_PATTERNS.${processParams.off} \
		${processParams?.c ? "-c " + processParams.c  : ""} \
		${processParams?.l ? "-l " + processParams.l  : ""} \
		${processParams?.estimate_l ? "--estimate-l " + processParams.estimate_l  : ""} \
		${processParams?.adjust_l ? "--adjust-l " + processParams.adjust_l  : ""} \
		--ncpu ${task.cpus}
      	"""
}


process SPATIALDE__ADD_SPATIAL_PATTERNS {

    container params.tools.spatialde.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId), path(f2)

    output:
        tuple val(sampleId), path("${sampleId}.SPATIALDE__ADDED_SPATIAL_PATTERNS.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.spatialde)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/add_patterns_as_metric.py \
		${f} \
		${f2} \
      		${sampleId}.SPATIALDE__ADDED_SPATIAL_PATTERNS.${processParams.off}
      	"""
}
