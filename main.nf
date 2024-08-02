#!/usr/bin/env nextflow
include { BasecallingAndDemux } from './subworkflows/basecalling_demux.nf'
include { QualityCheck }        from './subworkflows/quality_check.nf'
include { GenerateReports }     from './subworkflows/reports.nf'
include { CollectVersions }     from './subworkflows/versions.nf'

include { pathCheck } from './lib/groovy/utils.gvy'


// check and prepare input channels
data_dir = pathCheck(params.data_dir, isDirectory = true)
multiqc_config = pathCheck("${workflow.projectDir}/tool_conf/multiqc_config.yaml")

if (params.skip_demultiplexing) {
  sample_names = channel.fromList([])
} else {
  pathCheck(params.sample_data)
  channel
    .fromPath(params.sample_data)
    .splitCsv(header: true)
    .map { row -> [row.barcode, row.sample] }
    .set { sample_names }
}

workflow {
  BasecallingAndDemux(sample_names, data_dir)

  QualityCheck(
    BasecallingAndDemux.out.sequences
  )

  GenerateReports(
    QualityCheck.out.software_reports,
    BasecallingAndDemux.out.sequencing_summary,
    data_dir,
    sample_names.map { it[0] }.collect(),
    multiqc_config
  )
}