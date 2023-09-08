nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*
import org.yaml.snakeyaml.Yaml

include {
	isParamNull;
} from '../../utils/processes/utils.nf'

/* 
 * STATIC VERSION GENERATE REPORT
 * 
 * General reporting function: 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
 */
process SPATIAL__SQUIDPY__GENERATE_REPORT {

  	container params.tools.squidpy.container
  	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
    label 'compute_resources__report'

	input:
		file ipynb
		tuple val(sampleId), path(adata)
		val(reportTitle)

	output:
		tuple val(sampleId), path("${sampleId}.${reportTitle}.ipynb")

	script:
		def reportParams = new Yaml().dump(annotations_to_plot: params.tools.squidpy.report.annotations_to_plot)
		"""
		papermill ${ipynb} \
		    --report-mode \
			${sampleId}.${reportTitle}.ipynb \
			-p FILE $adata \
			-y "${reportParams}" \
			-p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
			-p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
		"""

}

process SPATIAL__SQUIDPY__REPORT_TO_HTML {

	container params.tools.squidpy.container
	publishDir "${params.global.outdir}/notebooks/intermediate", mode: 'link', overwrite: true
	// copy final "merged_report" to notbooks root:
	publishDir "${params.global.outdir}/notebooks", pattern: '*merged_report*', mode: 'link', overwrite: true
    label 'compute_resources__report'

	input:
		tuple val(sampleId), path(ipynb)

	output:
		file("*.html")

	script:
		"""
		jupyter nbconvert ${ipynb} --to html
		"""

}
