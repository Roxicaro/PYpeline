#!/usr/bin/env nextflow

process SAMTOOLS_INDEXING {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    conda "envs/bwa.yaml"
    container "community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509"

    input:
    tuple val(sample_id), path(sorted)

    output:
    path "${sample_id}.bam.bai"

    script:
    """
    samtools index ${sorted}
    """
}