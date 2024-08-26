include { slugify } from '../lib/groovy/utils.gvy'


workflow QualityCheck {
  take:
    sequences   // channel [name, fastq]

  main:
    sequences
      | (fastQC & nanoPlot & nanoq)
      | mix
      | map { it[1] }
      | collect
      | set { reports }

  emit:
    software_reports = reports
}


process fastQC {
  label 'fastqc'
  tag { name }
  publishDir "${params.output_dir}/qc/fastqc", mode: 'copy'
  cpus { 4 * task.attempt }
  memory { 8.GB * task.attempt }
  errorStrategy 'retry'
  maxRetries 3

  when:
  'fastqc' in params.qc_tools

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
  tag { name }
  publishDir "${params.output_dir}/qc/nanoplot", mode: 'copy'
  cpus 4

  when:
  'nanoplot' in params.qc_tools

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


process nanoq {
  label 'nanoq'
  tag { name }
  publishDir "${params.output_dir}/qc/nanoq", mode: 'copy'
  cpus 1

  when:
  'nanoq' in params.qc_tools

  input:
  tuple val(name), path(reads)
  
  output:
  tuple val(name), path("${name}.txt")
  
  script:
  """
  nanoq -svv -i ${reads} > ${name}.txt
  """
}
