include { sanitizeFilename } from '../lib/groovy/utils.gvy'


workflow QualityCheck {
  take:
    sequences           // channel [name, fastq]
    sequencing_summary  // sequencing summary file

  main:
    pycoQC(sequencing_summary)

    sequences
      | (fastQC & nanoPlot)
      | mix
      | map { it[1] }
      | collect
      | set { reports }

  emit:
    software_reports = reports
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
  file_opt = reads.name.endsWith('.ubam') ? '--ubam' : '--fastq'
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
