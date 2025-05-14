#!/usr/bin/env nextflow
include { validateParameters ; samplesheetToList } from 'plugin/nf-schema'
include { BasecallingAndDemux } from './subworkflows/basecalling_demux.nf'
include { QualityCheck        } from './subworkflows/quality_check.nf'
include { GenerateReports     } from './subworkflows/reports.nf'


// validate and prepare input channels
workflow {
  validateParameters()

  data_dir = file(params.data_dir, type: 'dir', checkIfExists: true)
  multiqc_config = file("${workflow.projectDir}/tool_conf/multiqc_config.yaml", checkIfExists: true)

  if (params.sample_data) {
    samples = channel
      .fromList(samplesheetToList(params.sample_data, "assets/samples_data_schema.json"))
      .map { row -> [[id: row[1], barcode: row[0]]] }
  }
  else {
    samples = channel.empty()
  }


  BasecallingAndDemux(samples, data_dir)

  QualityCheck(BasecallingAndDemux.out.sequences)

  GenerateReports(
    QualityCheck.out.software_reports,
    BasecallingAndDemux.out.sequencing_summary,
    samples.map { row -> row[0].barcode }.collect().ifEmpty { [] },
    data_dir,
    multiqc_config,
  )
}
