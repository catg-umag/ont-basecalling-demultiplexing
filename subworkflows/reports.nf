include { sanitizeFilename } from '../lib/groovy/utils.gvy'


workflow GenerateReports {
  take:
    software_reports    // channel [file] | SW reports
    software_versions   // channel [file] | SW versions file
    model_versions      // channel [file] | SW model versions file
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