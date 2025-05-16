workflow BasecallingAndDemux {
    take:
    sample_names // channel [barcode, sample]
    data_dir     // directory containing POD5 files

    main:
    basecalling(data_dir)
    getDoradoModel(basecalling.out.reads.map { it[1] })
    qscoreFiltering(basecalling.out.reads)

    ch_sequences_for_qc = qscoreFiltering.out.reads_pass.mix(qscoreFiltering.out.reads_fail)

    if (params.sample_data) {
        /* Demultiplexing */
        demultiplexing(qscoreFiltering.out.reads_pass)

        ch_classified_reads = demultiplexing.out.classified
            .flatMap { exp ->
                exp.collect { file ->
                    [(file.name =~ /barcode[0-9]+/)[0], file]
                }
            }
            .join(sample_names.map { sample -> [sample[0].barcode] + sample })
            .map { [it[2], it[1]] }
            .map { meta, file -> [meta + [step: 'demux', classified: true], file] }

        ch_sequences_for_qc = ch_sequences_for_qc
            .mix(ch_classified_reads)
            .mix(
                demultiplexing.out.unclassified.map { file ->
                    [[id: 'unclassified', step: 'demux', classified: false], file]
                }
            )

        mergeDoradoSummaries(
            basecalling.out.sequencing_summary,
            demultiplexing.out.summary,
        )
        ch_sequencing_summary = mergeDoradoSummaries.out.summary
    }
    else {
        ch_sequencing_summary = basecalling.out.sequencing_summary
    }

    // If FASTQ output is disabled, rename the BAM files
    if (!params.fastq_output) {
        renameUbam(ch_sequences_for_qc)
    }

    // QC tools only support FASTQ files, so convertion is always neccessary
    bamToFastq(ch_sequences_for_qc)

    emit:
    sequences          = bamToFastq.out.fastq
    sequencing_summary = ch_sequencing_summary
}


process basecalling {
    label 'dorado'
    clusterOptions "--gres=gpu:${params.dorado_basecalling_gpus}"
    cpus { 4 * params.dorado_basecalling_gpus }
    memory "${16 * params.dorado_basecalling_gpus}G"

    input:
    path data_dir

    output:
    tuple val { [id: 'basecalled'] }, path('basecalled.ubam'), emit: reads
    path 'sequencing_summary.txt', emit: sequencing_summary
    tuple val('dorado'), eval("dorado --version 2>&1 | tr -d '\n'"), topic: versions

    script:
    extra_args = task.ext.args ?: ''
    """
    dorado basecaller \\
        --recursive \\
        --device 'cuda:all' \\
        --trim adapters \\
        ${extra_args} \\
        ${params.dorado_basecalling_model} \\
        ${data_dir} \\
    > basecalled.ubam

    dorado summary basecalled.ubam > sequencing_summary.txt
    """
}


process qscoreFiltering {
    label 'samtools'
    cpus 4

    input:
    tuple val(meta), path(ubam)

    output:
    tuple val { [id: 'basecalled_pass', step: 'basecall'] }, path('basecalled_pass.ubam'), emit: reads_pass
    tuple val { [id: 'basecalled_fail', step: 'basecall'] }, path('basecalled_fail.ubam'), emit: reads_fail
    tuple val('samtools'), eval("samtools --version | head -n 1 | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    samtools view \\
        -e '[qs] >= ${params.qscore_filter}' ${ubam} \\
        --output basecalled_pass.ubam \\
        --unoutput basecalled_fail.ubam \\
        --bam \\
        --threads ${task.cpus} 
    """
}


process demultiplexing {
    label 'dorado'
    cpus 8

    input:
    tuple val(meta), path(basecalled_reads)

    output:
    path 'demultiplexed/*barcode*', emit: classified
    path 'demultiplexed/*unclassified*', emit: unclassified
    path 'demultiplexed/barcoding_summary.txt', emit: summary
    tuple val('dorado'), eval("dorado --version 2>&1 | tr -d '\n'"), topic: versions

    script:
    both_ends = params.dorado_demux_both_ends ? '--barcode-both-ends' : ''
    extra_args = task.ext.args ?: ''
    """
    dorado demux \\
        --output-dir demultiplexed/ \\
        --kit-name ${params.dorado_demux_kit} \\
        --emit-summary \\
        --threads ${task.cpus} \\
        ${both_ends} \\
        ${extra_args} \\
        ${basecalled_reads}
    """
}


process renameUbam {
    label 'linux'
    tag { meta.id }

    input:
    tuple val(meta), path(bam, stageAs: 'reads.orig.ubam')

    output:
    tuple val(meta), path("${meta.id}.ubam")

    script:
    """
    cp ${bam} ${meta.id}.ubam
    """
}


process bamToFastq {
    label 'samtools'
    tag { meta.id }
    cpus 4

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq
    tuple val('samtools'), eval("samtools --version | head -n 1 | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    samtools fastq -T '*' -@ ${task.cpus} -0 ${meta.id}.fastq.gz ${bam}
    """
}


process getDoradoModel {
    label 'samtools'

    input:
    path bam

    output:
    path 'dorado_model.tsv', topic: model_versions
    tuple val('samtools'), eval("samtools --version | head -n 1 | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    model_version=\$(samtools view -H ${bam} | grep -Po '(?<=basecall_model=)([^ ]+)' | uniq)
    echo "Software\tModel\tVersion" > dorado_model.tsv
    echo "Dorado\tBasecalling\t\${model_version}" >> dorado_model.tsv 
    """
}


process mergeDoradoSummaries {
    label 'csvtk'

    input:
    path sequencing_summary, stageAs: 'input/sequencing_summary.txt'
    path barcoding_summary, stageAs: 'input/barcoding_summary.txt'

    output:
    path 'sequencing_summary.txt', emit: summary
    tuple val('csvtk'), eval("csvtk version | grep -Eo '[0-9.]+' | tr -d '\n'"), topic: versions

    script:
    """
    # Drop barcode column from sequencing_summary,
    # add barcode column from barcoding_summary
    # and fix empty barcode values
    csvtk cut -t \\
        --fields -barcode \\
        ${sequencing_summary} \\
    | csvtk join -t \\
        --left-join \\
        --fields 'filename,read_id' \\
        - ${barcoding_summary} \\
    | csvtk replace -t \\
        --fields barcode \\
        --pattern '^\$' \\
        --replacement 'unclassified' \\
        --out-file sequencing_summary.txt
    """
}
