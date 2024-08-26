#!/usr/bin/env nextflow
include { BasecallingAndDemux } from './subworkflows/basecalling_demux.nf'
include { QualityCheck }        from './subworkflows/quality_check.nf'
include { GenerateReports }     from './subworkflows/reports.nf'
include { CollectVersions }     from './subworkflows/versions.nf'


// check and prepare input channels
data_dir = file(params.data_dir, checkIfExists: true, type: 'dir')
multiqc_config = file("${workflow.projectDir}/tool_conf/multiqc_config.yaml", checkIfExists: true)

if (params.skip_demultiplexing) {
  samples = channel.fromList([])
} else {
  file(params.sample_data, checkIfExists: true)
  channel
    .fromPath(params.sample_data)
    .splitCsv(header: true)
    .map { row -> [row.barcode, row.sample] }
    .set { samples }
}

params.qc_tools = params.qc_tools.each { it.toLowerCase() }


workflow {
  BasecallingAndDemux(samples, data_dir)

  QualityCheck(BasecallingAndDemux.out.sequences)

  CollectVersions(BasecallingAndDemux.out.basecalled_ubam)

  GenerateReports(
    QualityCheck.out.software_reports,
    CollectVersions.out.software_versions,
    CollectVersions.out.model_versions,
    BasecallingAndDemux.out.sequencing_summary,
    samples.map { it[0] }.collect(),
    data_dir,
    multiqc_config
  )
}
