workflow BasecallingAndDemux {
  take:
    sample_names  // channel [barcode, sample]
    data_dir      // directory containing POD5 files

  main:
    basecalling(
      data_dir,
      downloadModel(params.dorado_basecalling_model).out
    )

    if (params.skip_demultiplexing && params.fastq_output) {
      bamToFastq(basecalling.out.bam_pass)
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
  publishDir "${params.output_dir}/bam/", \
    pattern: 'basecalled_pass.bam', \
    mode: 'copy', \
    enabled: params.skip_demultiplexing && !params.fastq_output
  clusterOptions = "--gres=gpu:${params.dorado_basecalling_gpus}"
  cpus params.dorado_basecalling_gpus
  
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
    ${params.guppy_basecalling_extra_config} \
    ${model} \
    ${data_dir} \
  > basecalled_pass.bam

  dorado summary basecalled_pass.bam > sequencing_summary.txt
  """
}


process demultiplexing {
  label 'dorado'
  cpus params.dorado_demux_cpus

  input:
  path(basecalled_reads)

  output:
  path('demultiplexed/*barcode*')                , emit: classified
  path('demultiplexed/unclassified*', arity: '1'), emit: unclassified

  script:
  both_ends = params.dorado_demux_both_ends ? '--barcode-both-ends' : ''
  emit_fastq = params.fastq_output ? '--emit-fastq' : ''
  """
  dorado demux \
    --output-dir demultiplexed/ \
    --kit-name "${params.dorado_demux_kit}" \
    ${both_ends} \
    ${emit-fastq} \
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
  tuple val(name), path("${name}.{fastq.gz,bam}", arity: '1')
  
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
  publishDir "${params.output_dir}/fastq/", mode: 'copy'
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