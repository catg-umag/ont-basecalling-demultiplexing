include { sanitizeFilename } from '../lib/groovy/utils.gvy'


workflow GenerateReports {
  take:
    software_reports    // channel [file] | SW reports
    software_versions   // channel [file] | SW versions file
    model_versions      // channel [file] | SW model versions file
    sequencing_summary  // channel [file] | Sequencing summary file
    multiqc_config      // multiqc config file

  main:
    multiQC(software_reports, software_versions, model_versions, multiqc_config)
}


process multiQC {
  label 'multiqc'
  publishDir "${params.output_dir}/multiqc", mode: 'copy'
  
  input:
  path(reports, stageAs: 'reports/*')
  path('reports/versions/*_mqc_versions.yaml')
  path('reports/model_versions/model_versions_*.tsv')
  path('multiqc_config.yaml')

  output:
  tuple path('*multiqc_data'), path('*multiqc*.html')

  script:
  if (params.experiment_name) {
    filename = sanitizeFilename("${params.experiment_name}_multiqc")
    title_opts = "--title '${params.experiment_name} Report' --filename ${filename}"
  } else {
    title_opts = ''
  }
  """
  multiqc ${title_opts} reports
  """
}


process toulligQC {
  label 'toulligqc'
  publishDir "${params.output_dir}/reports/toulligqc", mode: 'copy'
  memory 4.GB

  input:
  path(data_dir)
  path(sequencing_summary)
  val(barcode_list)
  
  output:
  path(report_filename)
  tuple val('ToulligQC'), eval('toulligqc --version'), topic: versions
  
  script:
  report_filename = "${sanitizeFilename(params.experiment_name)}_toulligqc.html"
  barcodes_opt = barcode_list
    ? "--barcoding --barcodes ${barcode_list.join(',')}"
    : ''
  """
  toulligqc \
    --sequencing-summary-source ${sequencing_summary} \
    --pod5-source ${data_dir} \
    --html-report-path ${report_filename} \
    --report-name ${params.experiment_name} \
    --qscore-threshold ${params.qscore_filter} \
    ${barcodes_opt}
  """
}


process pycoQC {
  label 'pycoqc'
  publishDir "${params.output_dir}/qc/pycoqc", mode: 'copy'
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 5

  input:
  path(sequencing_summary)
  
  output:
  path('pycoQC_report.html')
  
  script:
  title_opt = params.experiment_name
    ? "--report_title '${params.experiment_name} Sequencing Report'"
    : ''
  """
  pycoQC \
    -f ${sequencing_summary} \
    ${title_opt} \
    -o pycoQC_report.html
  """
}