#!/usr/bin/env nextflow

process SAMTOOLS_INDEXING {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    conda "envs/bwa.yaml"

    input:
    tuple val(sample_id), path(sorted)

    output:
    path "${sample_id}.bam.bai"

    script:
    """
    samtools index ${sorted}
    """
}