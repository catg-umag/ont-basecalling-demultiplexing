include { slugify } from '../lib/nf/utils.nf'


workflow QualityCheck {
    take:
    sequences // channel [[id: 'id', ...], fastq]

    main:
    ch_reports = channel.empty()

    if ('fastqc' in params.qc_tools) {
        fastQC(sequences)
        ch_reports = ch_reports.mix(fastQC.out.reports)
    }
    if ('nanoplot' in params.qc_tools) {
        nanoPlot(sequences)
        ch_reports = ch_reports.mix(nanoPlot.out.reports)
    }

    if ('nanoq' in params.qc_tools) {
        nanoq(sequences)
        ch_reports = ch_reports.mix(nanoq.out.reports)
    }

    emit:
    software_reports = ch_reports.map { it[1] }
}


process fastQC {
    label 'fastqc'
    tag { meta.id }
    cpus { task.attempt == 1 ? 4 : 8 }
    memory { 4.GB * Math.pow(2, task.attempt - 1) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("fastqc_${meta.id}"), emit: reports
    tuple val('FastQC'), eval("fastqc --version | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    mkdir fastqc_${meta.id}
    fastqc \
        ${reads} \
        -o fastqc_${meta.id} \
        -t ${task.cpus} \
        --memory ${task.memory.toMega() / task.cpus}
    """
}


process nanoPlot {
    label 'nanoplot'
    tag { meta.id }
    cpus 4
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("nanoplot_${meta.id}"), emit: reports
    tuple val('NanoPlot'), eval("NanoPlot --version | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    NanoPlot \
        --fastq ${reads} \
        --outdir nanoplot_${meta.id} \
        --prefix ${meta.id}_ \
        --threads ${task.cpus}
    """
}


process nanoq {
    label 'nanoq'
    tag { meta.id }
    cpus 1
    memory 2.GB

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: reports
    tuple val('nanoq'), eval("nanoq --version | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    nanoq -svv -i ${reads} > ${meta.id}.txt
    """
}
