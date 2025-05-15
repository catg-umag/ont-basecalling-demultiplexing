#!/usr/bin/env nextflow
include { validateParameters ; samplesheetToList } from 'plugin/nf-schema'
include { BasecallingAndDemux } from './subworkflows/basecalling_demux.nf'
include { QualityCheck        } from './subworkflows/quality_check.nf'
include { GenerateReports     } from './subworkflows/reports.nf'
include { SequenceCollection  } from './subworkflows/sequence_collection.nf'


workflow {
  // validate and prepare input channels
  validateParameters()

  if (params.basecalled_input) {
    if ('toulligqc' in params.qc_tools) {
      log.warn("ToulligQC is not supported for basecalled input.")
    }
    if ('pycoqc' in params.qc_tools) {
      log.warn("pycoQC is not supported for basecalled input.")
    }
    if (!params.fastq_output) {
      log.warn("BAM output not supported with basecalled input.")
    }
  }

  data_dir = file(params.data_dir, type: 'dir', checkIfExists: true)
  multiqc_config = file("${workflow.projectDir}/tool_conf/multiqc_config.yaml", checkIfExists: true)

  if (params.sample_data) {
    samples = Channel.fromList(samplesheetToList(params.sample_data, "assets/samples_data_schema.json"))
      .map { row -> [[id: row[1], barcode: row[0]]] }
  }
  else {
    samples = Channel.empty()
  }


  if (params.basecalled_input) {
    SequenceCollection(samples, data_dir)
    sequences = SequenceCollection.out.sequences
    sequencing_summary = Channel.empty()
  }
  else {
    BasecallingAndDemux(samples, data_dir)
    sequences = BasecallingAndDemux.out.sequences
    sequencing_summary = BasecallingAndDemux.out.sequencing_summary
  }

  QualityCheck(sequences)

  GenerateReports(
    QualityCheck.out.software_reports,
    sequencing_summary,
    samples.map { row -> row[0].barcode }.collect().ifEmpty { [] },
    data_dir,
    multiqc_config,
  )
}
