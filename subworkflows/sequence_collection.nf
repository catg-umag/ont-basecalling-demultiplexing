workflow SequenceCollection {
  take:
  sample_names // channel [barcode, sample]
  data_dir     // directory containing POD5 files

  main:
  if (params.sample_data) {
    ch_barcode_files = Channel.from(data_dir.listFiles())
      .map { dir -> [dir.name, dir] }
      .join(sample_names.map { sample -> [sample[0].barcode] + sample }, remainder: true)
      .branch { _barcode, dir, meta ->
        classified: meta
        return [meta, dir]
        unclassified: true
        return [[id: 'unclassified'], dir]
      }
    ch_sequences_to_collect = ch_barcode_files.classified.mix(
      ch_barcode_files.unclassified.groupTuple()
    )
  }
  else {
    ch_sequences_to_collect = Channel.from([[['id': 'basecalled'], data_dir]])
  }

  concatFastq(ch_sequences_to_collect)

  emit:
  sequences = concatFastq.out.sequences
}


process concatFastq {
  label 'linux'
  tag { meta.id }

  input:
  tuple val(meta), path("input/*")

  output:
  tuple val(meta), path("${meta.id}.fastq.gz"), emit: sequences

  script:
  """
  files=\$(find -L input -type f -name "*.fastq.gz")
  cat \${files} > ${meta.id}.fastq.gz
  """
}
