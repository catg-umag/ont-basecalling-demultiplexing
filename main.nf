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
  sample_names = channel
    .fromPath(params.sample_data)
    .splitCsv(header: true)
    .map { row -> [row.barcode, row.sample] }
}

workflow {
  BasecallingAndDemux(sample_names, data_dir)

  QualityCheck(
    BasecallingAndDemux.out.sequences,
    BasecallingAndDemux.out.sequencing_summary
  )

  CollectVersions()

  GenerateReports(
    QualityCheck.out.software_reports,
    CollectVersions.out.software_versions,
    CollectVersions.out.model_versions,
    multiqc_config
  )
}