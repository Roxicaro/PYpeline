#!/usr/bin/env nextflow

process SAMTOOLS_SORTING {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    conda "envs/bwa.yaml"

    input:
    tuple val(sample_id), path(mapped)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools sort -o "${sample_id}.bam" ${mapped}
    """
}