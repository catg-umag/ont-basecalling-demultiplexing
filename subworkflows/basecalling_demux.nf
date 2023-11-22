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
      demultiplexing(basecalling.out.fastq_pass)

      demultiplexing.out.classified
        | flatMap { it.collect { x -> [x.simpleName, x] } }
        | join(sample_names)
        | map { [it[2], it[1]] }
        | mix(demultiplexing.out.unclassified.map { ['unclassified', it] })
        | mix(sequences_to_merge)
        | mergeSequences
        | filter { !['unclassified', 'pass', 'fail'].contains(it[0]) }
        | compressSequences
    }

  emit:
    sequences = mergeSequences.out
    sequencing_summary = basecalling.out.sequencing_summary
    barcoding_summary = params.skip_demultiplexing 
      ? file('NO_FILE')
      : demultiplexing.out.barcoding_summary
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
    enabled: { params.skip_demultiplexing && !params.fastq_output }
  clusterOptions = "--gres=gpu:${params.dorado_basecalling_gpus}"
  cpus params.dorado_basecalling_gpus
  
  input:
  path(data_dir)
  path(basecalling_model)

  output:
  tuple val('basecalled_pass'), path('basecalled_pass.bam'), emit: bam_pass
  path('sequencing_summary.txt'), emit: sequencing_summary

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
  publishDir "${params.output_dir}/sequencing_info/", \
    pattern: 'demultiplexed/barcoding_summary.txt', \
    saveAs: { 'barcoding_summary.txt' }, \
    mode: 'copy'
  cpus params.guppy_barcoding_cpus

  input:
  path(fastq_dir)

  output:
  path('demultiplexed/barcode*'), emit: classified
  path('demultiplexed/unclassified'), emit: unclassified
  path('demultiplexed/barcoding_summary.txt'), emit: barcoding_summary

  script:
  both_ends = params.guppy_barcoding_both_ends ? '--require_barcodes_both_ends' : ''
  """
  guppy_barcoder \
    --input_path ${fastq_dir} \
    --save_path demultiplexed/ \
    --recursive \
    --barcode_kits "${params.guppy_barcoding_kits}" \
    ${both_ends} \
    --enable_trim_barcodes \
    --detect_adapter \
    --detect_barcodes \
    --worker_threads ${task.cpus} \
    ${params.guppy_barcoding_extra_config}
  """
}


process mergeSequences {
  label 'linux'
  publishDir "${params.output_dir}/basecalled/", \
    pattern: '*.fastq', \
    mode: 'copy', \
    enabled: params.skip_demultiplexing
  tag "${name}"

  input:
  tuple val(name), path(fastq_dir)

  output:
  tuple val(name), path("${name}.fastq")

  script:
  """
  cat ${fastq_dir}/*.fastq > ${name}.fastq
  """
}


process compressSequences {
  label 'pigz'
  tag "${name}"
  publishDir "${params.output_dir}/fastq/", mode: 'copy'
  cpus 4

  input:
  tuple val(name), path(fastq)
  
  output:
  tuple val(name), path("${name}.fastq.gz")
  
  script:
  """
  pigz -p ${task.cpus} -c ${fastq} > ${name}.fastq.gz
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