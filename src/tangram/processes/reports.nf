nextflow.enable.dsl=2

import static groovy.json.JsonOutput.*
import org.yaml.snakeyaml.Yaml

include {
        isParamNull;
} from '../../utils/processes/utils.nf'

/* 
 * takes a template ipynb and adata as input,
 * outputs ipynb named by the value in ${reportTitle}
 */
process TANGRAM__GENERATE_DUAL_INPUT_REPORT {
        def reportParams = new Yaml().dump(annotations_to_plot: params.tools.scanpy.report.annotations_to_plot)
	container params.tools.tangram.container
	
        publishDir "${params.global.outdir}/notebooks", mode: 'link', overwrite: true
    	label 'compute_resources__report'

        input:
                file(ipynb)
                tuple \
			val(sampleId), \
                        file(data1), \
                        file(data2), \
                        val(stashedParams)
                val(reportTitle)
                val(isParameterExplorationModeOn)

        output:
        tuple \
                        val(sampleId), \
                        file("${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + "." : ''}ipynb"), \
                        val(stashedParams)

        script:
                if(!isParamNull(stashedParams))
                        uuid = stashedParams.findAll { it != 'NULL' }.join('_')
                """
                papermill ${ipynb} \
                    --report-mode \
                        ${sampleId}.${reportTitle}.${isParameterExplorationModeOn ? uuid + "." : ''}ipynb \
                        -p FILE1 $data1 -p FILE2 $data2 \
                        -y "${reportParams}" \
                        -p WORKFLOW_MANIFEST '${params.misc.manifestAsJSON}' \
                        -p WORKFLOW_PARAMETERS '${params.misc.paramsAsJSON}'
                """		
}


process TANGRAM__REPORT_TO_HTML {

	container params.tools.scanpy.container
	publishDir "${params.global.outdir}/notebooks", mode: 'link', overwrite: true
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
