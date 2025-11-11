#!/usr/bin/env nextflow

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    conda "envs/fastqc.yaml"
    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"

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