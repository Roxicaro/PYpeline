#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

workflow {
    Channel.fromFilePairs(params.reads)
        .set { raw_reads }

    FASTQC(raw_reads)
}