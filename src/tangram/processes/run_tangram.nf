nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/tangram/bin/" : ""

process TANGRAM__MAP_CELLTYPES {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__gpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId2), path(ref)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__CELLTYPES.${processParams.off}")
	tuple val(sampleId), path("${sampleId}.TANGRAM__MAPPING.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/run_tangram.py \
		-r  ${ref} \
		--output-mapping ${sampleId}.TANGRAM__MAPPING.${processParams.off} \
		${processParams?.device ? "-d " + processParams.device  : ""} \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${processParams?.gene_selection_method ? "-m " + processParams.gene_selection_method  : ""} \
		${processParams?.number_genes ? "-n " + processParams.number_genes  : ""} \
		${processParams?.qvalue ? "-q " + processParams.qvalue  : ""} \
		${processParams?.max_sparsity ? "-s " + processParams.max_sparsity  : ""} \
		${processParams?.mapping_mode ? "--mode " + processParams.mapping_mode  : ""} \
		${processParams?.exp_scrnaseq ? "--exp-ref"  : ""} \
		${f} \
      		${sampleId}.TANGRAM__CELLTYPES.${processParams.off}
      	"""
}


process TANGRAM__PREPARE_SCRNA {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__SCRNA.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/prepare_scRNA_tangram.py \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${processParams?.gene_selection_method ? "-m " + processParams.gene_selection_method  : ""} \
		${processParams?.rank_gene_method ? "--rank-gene-method " + processParams.rank_gene_method  : ""} \
		${processParams?.normalize_scRNAseq ? "--normalize " + processParams.normalize_scRNAseq  : ""} \
		${processParams?.max_cells ? "--max-cells-cluster " + processParams.max_cells  : ""} \
		${params.global.seed ? "--seed " + params.global.seed : ""} \
		${f} \
      		${sampleId}.TANGRAM__SCRNA.${processParams.off}
      	"""
}


process TANGRAM__PREPARE_SPATIAL {

    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId), path(filtered)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__SPATIAL.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/prepare_spatial_tangram.py \
		${f} \
		${filtered} \
		${processParams?.normalize_spatial ? "--normalize " + processParams.normalize_spatial  : ""} \
      		${sampleId}.TANGRAM__SPATIAL.${processParams.off}
      	"""
}

process TANGRAM__PREPARE_SPATIAL_SIMPLE {

    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__SPATIAL.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/prepare_spatial_tangram_simple.py \
		${f} \
		${processParams?.normalize_spatial ? "--normalize " + processParams.normalize_spatial  : ""} \
      		${sampleId}.TANGRAM__SPATIAL.${processParams.off}
      	"""
}


process TANGRAM__COMPUTE_MAPPING {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__gpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId2), path(ref)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__MAPPING.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/map_dataset_tangram.py \
		-r  ${ref} \
		${processParams?.device ? "-d " + processParams.device  : ""} \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${processParams?.gene_selection_method ? "-m " + processParams.gene_selection_method  : ""} \
		${processParams?.number_genes ? "-n " + processParams.number_genes  : ""} \
		${processParams?.qvalue ? "-q " + processParams.qvalue  : ""} \
		${processParams?.max_sparsity ? "-s " + processParams.max_sparsity  : ""} \
		${processParams?.mapping_mode ? "--mode " + processParams.mapping_mode  : ""} \
		${processParams?.gene_list ? "--list_genes " + processParams.gene_list  : ""} \
		${params.global.seed ? "--seed " + params.global.seed : ""} \
		${f} \
      		${sampleId}.TANGRAM__MAPPING.${processParams.off}
      	"""
}


process TANGRAM__PROJECT_CELLTYPES {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId2), path(ref)
	tuple val(sampleId), path(mapping)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__CELLTYPES.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/project_celltype_tangram.py \
		-r  ${ref} \
		--mapping ${mapping} \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${processParams?.normalize_celltype_scores ? "--normalize_map_scores " + processParams.normalize_celltype_scores  : ""} \
		${processParams?.quantile_mapping_score ? "--quantile-score " + processParams.quantile_mapping_score  : ""} \
		${f} \
      		${sampleId}.TANGRAM__CELLTYPES.${processParams.off}
      	"""
}


process TANGRAM__PROJECT_EXPRESSION {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)
	tuple val(sampleId2), path(ref)
	tuple val(sampleId), path(mapping)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__GEX.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/project_expression_tangram.py \
		-r  ${ref} \
		--mapping ${mapping} \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${f} \
      		${sampleId}.TANGRAM__GEX.${processParams.off}
      	"""
}
