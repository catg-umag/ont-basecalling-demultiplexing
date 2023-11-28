workflow BasecallingAndDemux {
  take:
    sample_names  // channel [barcode, sample]
    data_dir      // directory containing POD5 files

  main:
    basecalling(
      data_dir,
      downloadModel(params.dorado_basecalling_model)
    )

    if (params.skip_demultiplexing) {
      if (params.fastq_output) {
        bamToFastq(basecalling.out.reads_pass)
          | set { sequences_to_postprocess }
      } else {
        sequences_to_postprocess = basecalling.out.reads_pass
      }
    } else {
      demultiplexing(basecalling.out.reads_pass)

      demultiplexing.out.classified
        | flatMap { it.collect { x -> [(x.name =~ /barcode[0-9]+/)[0], x] } }
        | join(sample_names)
        | map { [it[2], it[1]] }
        | gatherSequences

      gatherSequences.out
        | mix(demultiplexing.out.unclassified.map { ['unclassified', it] })
        | set { sequences_to_postprocess }
    }

  emit:
    sequences = sequences_to_postprocess
    sequencing_summary = basecalling.out.sequencing_summary
}


process downloadModel {
  label 'dorado'

  input:
  val(model_name)
  
  output:
  path(model_name)
  
  script:
  """
  dorado download --model ${model_name}
  """
}


process basecalling {
  label 'dorado'
  publishDir "${params.output_dir}/sequencing_info/", \
    pattern: 'sequencing_summary.txt', \
    mode: 'copy'
  publishDir "${params.output_dir}/basecalled/", \
    pattern: 'basecalled_pass.bam', \
    mode: 'copy', \
    enabled: params.skip_demultiplexing && !params.fastq_output
  clusterOptions = "--gres=gpu:${params.dorado_basecalling_gpus}"
  cpus { 4 * params.dorado_basecalling_gpus }
  memory "${16 * params.dorado_basecalling_gpus}G"
  
  input:
  path(data_dir)
  path(basecalling_model)

  output:
  tuple val('basecalled_pass'), path('basecalled_pass.bam'), emit: reads_pass
  path('sequencing_summary.txt')                           , emit: sequencing_summary

  script:
  """
  dorado basecaller \
    --min-qscore ${params.dorado_basecalling_min_qscore} \
    --recursive \
    --device 'cuda:all' \
    ${params.dorado_basecalling_extra_config} \
    ${basecalling_model} \
    ${data_dir} \
  > basecalled_pass.bam

  dorado summary basecalled_pass.bam > sequencing_summary.txt
  """
}


process demultiplexing {
  label 'dorado'
  cpus params.dorado_demux_cpus

  input:
  tuple val(name), path(basecalled_reads)

  output:
  path('demultiplexed/*barcode*')    , emit: classified
  path('demultiplexed/unclassified*'), emit: unclassified

  script:
  both_ends = params.dorado_demux_both_ends ? '--barcode-both-ends' : ''
  emit_fastq = params.fastq_output ? '--emit-fastq' : ''
  """
  dorado demux \
    --output-dir demultiplexed/ \
    --kit-name "${params.dorado_demux_kit}" \
    ${both_ends} \
    ${emit_fastq} \
    --threads ${task.cpus} \
    ${params.dorado_demux_extra_config} \
    ${basecalled_reads}
  """
}


process gatherSequences {
  label 'pigz'
  tag "${name}"
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
  tag "${name}"
  publishDir "${params.output_dir}/basecalled/", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(bam)

  output:
  tuple val(name), path("${name}.fastq.gz")

  script:
  """
  samtools fastq -@ ${task.cpus} -0 ${name}.fastq.gz ${bam}
  """
}