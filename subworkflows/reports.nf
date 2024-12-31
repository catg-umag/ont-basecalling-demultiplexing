include { slugify } from '../lib/nf/utils.nf'


workflow GenerateReports {
  take:
  software_reports   // channel [file] | SW reports
  sequencing_summary // channel [file] | Sequencing summary file
  barcodes           // list [barcode]
  pod5_data          // directory containing POD5 files
  multiqc_config     // multiqc config file

  main:
  software_versions = channel
    .topic('versions')
    .unique()
    .collectFile(name: 'software_versions.yaml', newLine: true, sort: true) { tool, version ->
      "${tool}: \"${version}\""
    }
    .collect()
  model_versions = channel
    .topic('model_versions')
    .collect()

  if ('pycoqc' in params.qc_tools) {
    pycoQC(sequencing_summary)
  }
  if ('toulligqc' in params.qc_tools) {
    toulligQC(pod5_data, sequencing_summary, barcodes)
  }

  multiQC(software_reports.collect(), software_versions, model_versions, multiqc_config)
}


process toulligQC {
  label 'toulligqc'
  publishDir "${params.output_dir}/reports/toulligqc", mode: 'copy'
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
  maxRetries 3

  input:
  path data_dir
  path 'sequencing_summary.txt'
  val barcode_list

  output:
  path report_filename
  tuple val('ToulligQC'), eval('toulligqc --version'), topic: versions

  script:
  if (params.experiment_name) {
    report_filename = "${slugify(params.experiment_name)}_toulligqc.html"
    name_opt = "--report-name '${params.experiment_name}'"
  }
  else {
    report_filename = "toulligqc.html"
    name_opt = ''
  }
  barcodes_opt = params.sample_data && barcode_list
    ? "--barcoding --barcodes ${barcode_list.join(',')}"
    : ''
  """
  sed 's/${params.dorado_demux_kit}_//' sequencing_summary.txt > sequencing_summary.mod.txt

  toulligqc \
    --sequencing-summary-source sequencing_summary.mod.txt \
    --pod5-source ${data_dir} \
    --html-report-path ${report_filename} \
    --qscore-threshold ${params.qscore_filter} \
    ${name_opt} \
    ${barcodes_opt}
  """
}


process pycoQC {
  label 'pycoqc'
  publishDir "${params.output_dir}/reports/pycoqc", mode: 'copy'
  memory { 8.GB * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
  path sequencing_summary

  output:
  path 'pycoQC_report.html'
  tuple val('pycoQC'), eval("pycoQC --version | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

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


process multiQC {
  label 'multiqc'
  publishDir "${params.output_dir}/reports/multiqc", mode: 'copy'

  input:
  path reports, stageAs: 'reports/*'
  path 'reports/versions/ont_demux_*_mqc_versions.yaml'
  path 'reports/model_versions/model_versions_*.tsv'
  path 'multiqc_config.yaml'

  output:
  tuple path('*multiqc_data'), path('*multiqc*.html')

  script:
  if (params.experiment_name) {
    filename = slugify("${params.experiment_name}_multiqc")
    title_opts = "--title '${params.experiment_name} Report' --filename ${filename}"
  }
  else {
    title_opts = ''
  }
  """
  multiqc ${title_opts} reports
  """
}
