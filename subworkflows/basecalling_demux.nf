workflow BasecallingAndDemux {
  take:
    sample_names  // channel [barcode, sample]
    data_dir      // directory containing POD5 files

  main:
    basecalling(data_dir)
    qscoreFiltering(basecalling.out.reads)
      | mix
      | set { basecalled_reads_qscore_filtered }

    if (params.fastq_output) {
      basecalled_reads_qscore_filtered
        | bamToFastq
        | set { basecalled_reads }
    } else {
      basecalled_reads = basecalled_reads_qscore_filtered
    }

    if (params.sample_data) {
      demultiplexing(qscoreFiltering.out.reads_pass)

      demultiplexing.out.classified
        | flatMap { it.collect { x -> [(x.name =~ /barcode[0-9]+/)[0], x] } }
        | join(sample_names)
        | map { [it[2], it[1]] }
        | gatherSequences

      gatherSequences.out
        | mix(demultiplexing.out.unclassified.map { ['unclassified', it] })
        | mix(basecalled_reads)
        | set { sequences_to_postprocess }
    } else {
      sequences_to_postprocess = basecalled_reads
    }

  emit:
    basecalled_ubam = basecalling.out.reads.map { it[1] }
    sequences = sequences_to_postprocess
    sequencing_summary = basecalling.out.sequencing_summary
}


process basecalling {
  label 'dorado'
  publishDir "${params.output_dir}/sequencing_info/", \
    pattern: 'sequencing_summary.txt', \
    mode: 'copy'
  clusterOptions "--gres=gpu:${params.dorado_basecalling_gpus}"
  cpus { 4 * params.dorado_basecalling_gpus }
  memory "${16 * params.dorado_basecalling_gpus}G"
  
  input:
  path(data_dir)

  output:
  tuple val('basecalled'), path('basecalled.ubam'), emit: reads
  path('sequencing_summary.txt')                  , emit: sequencing_summary

  script:
  if (params.sample_data) {
    demux_opts = "--kit-name ${params.dorado_demux_kit}"
    demux_opts += params.dorado_demux_both_ends ? ' --barcode-both-ends' : ''
  } else {
    demux_opts = ''
  }
  extra_args = task.ext.args ?: ''
  """
  dorado basecaller \
    --recursive \
    --device 'cuda:all' \
    ${extra_args} \
    ${params.dorado_basecalling_model} \
    ${demux_opts} \
    ${data_dir} \
  > basecalled.ubam

  dorado summary basecalled.ubam > sequencing_summary.txt
  """
}


process qscoreFiltering {
  label 'samtools'
  publishDir "${params.output_dir}/basecalled/", \
    pattern: '*.bam', \
    mode: 'copy', \
    enabled: !params.sample_data && !params.fastq_output
  cpus 4

  input:
  tuple val(name), path(ubam)
  
  output:
  tuple val('reads_pass'), path('*.pass.ubam'), emit: reads_pass
  tuple val('reads_fail'), path('*.fail.ubam'), emit: reads_fail
  
  script:
  """
  samtools view \
    -e '[qs] >= ${params.qscore_filter}' ${ubam} \
    --output ${ubam.baseName}.pass.ubam \
    --unoutput ${ubam.baseName}.fail.ubam \
    --bam \
    --threads ${task.cpus} 
  """
}


process demultiplexing {
  label 'dorado'
  cpus 2

  input:
  tuple val(name), path(basecalled_reads)

  output:
  path('demultiplexed/*barcode*')    , emit: classified
  path('demultiplexed/*unclassified*'), emit: unclassified

  script:
  emit_fastq = params.fastq_output ? '--emit-fastq' : ''
  extra_args = task.ext.args ?: ''
  """
  dorado demux \
    --output-dir demultiplexed/ \
    --no-classify \
    ${emit_fastq} \
    --threads ${task.cpus} \
    ${extra_args} \
    ${basecalled_reads}
  """
}


process gatherSequences {
  label 'pigz'
  tag { name }
  publishDir "${params.output_dir}/demux/", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(reads)
  
  output:
  tuple val(name), path("${name}.{fastq.gz,bam}")
  
  script:
  if (params.fastq_output)
    """
    pigz -p ${task.cpus} -c ${reads} > ${name}.fastq.gz
    """
  else
    """
    cp ${reads} ${name}.bam
    """
}


process bamToFastq {
  label 'samtools'
  tag { name }
  publishDir "${params.output_dir}/basecalled/", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(bam)

  output:
  tuple val(name), path("${name}.fastq.gz")

  script:
  """
  samtools fastq -T '*' -@ ${task.cpus} -0 ${name}.fastq.gz ${bam}
  """
}
