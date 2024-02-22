workflow QualityCheck {
  take:
    sequences           // channel [name, fastq]
    sequencing_summary  // sequencing summary file
    multiqc_config      // multiqc config file

  main:
    pycoQC(sequencing_summary)

    sequences
      | (fastQC & nanoPlot)
      | mix
      | map { it[1] }
      | collect
      | multiMap {
          reports: it
          multiqc_config: multiqc_config
        }
      | multiQC
}


process fastQC {
  label 'fastqc'
  tag "${name}"
  publishDir "${params.output_dir}/qc/fastqc", mode: 'copy'
  cpus { 4 * task.attempt }
  memory { 8.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("fastqc_${name}")

  script:
  """
  mkdir fastqc_${name}
  fastqc \
    ${reads} \
    -o fastqc_${name} \
    -t ${task.cpus} --memory ${task.memory.toGiga()}GB
  """
}


process nanoPlot {
  label 'nanoplot'
  tag "${name}"
  publishDir "${params.output_dir}/qc/nanoplot", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), path("nanoplot_${name}")

  script:
  file_opt = reads.name.endsWith('.bam') ? '--ubam' : '--fastq'
  """
  NanoPlot \
    ${file_opt} ${reads} \
    --outdir nanoplot_${name} \
    --prefix ${name}_ \
    --threads ${task.cpus}
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


process multiQC {
  label 'multiqc'
  publishDir "${params.output_dir}/qc/multiqc", mode: 'copy'
  
  input:
  path(reports, stageAs: 'reports/*')
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