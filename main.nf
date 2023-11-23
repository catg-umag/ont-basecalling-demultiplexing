#!/usr/bin/env nextflow
include { addDefaultParamValues; pathCheck } from './lib/groovy/utils.gvy'

// load default parameters from YAML
addDefaultParamValues(params, "${workflow.projectDir}/params.default.yml")


include { BasecallingAndDemux } from './subworkflows/basecalling_demux.nf'
include { QualityCheck }        from './subworkflows/quality_check.nf'


// check and prepare input channels
data_dir = pathCheck(params.data_dir, isDirectory = true)
multiqc_config = pathCheck("${workflow.projectDir}/conf/multiqc_config.yaml")

if (params.skip_demultiplexing) {
  sample_names = channel.fromList([])
} else {
  pathCheck(params.sample_data)
  sample_names = channel
    .fromPath(params.sample_data)
    .splitCsv(header: true)
    .map { row -> [row.barcode, row.sample] }
}

workflow {
  BasecallingAndDemux(sample_names, data_dir)
  QualityCheck(
    BasecallingAndDemux.out.sequences,
    BasecallingAndDemux.out.sequencing_summary,
    multiqc_config
  )
}