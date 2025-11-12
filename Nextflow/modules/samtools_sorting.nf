#!/usr/bin/env nextflow

process SAMTOOLS_SORTING {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    conda "envs/bwa.yaml"
    container "community.wave.seqera.io/library/samtools_awscli:4fb4d9b382bd3990"
    //awscli	= 2.31.33   (conda-forge)
    //samtools  = 1.22.1    (bioconda)

    input:
    tuple val(sample_id), path(mapped)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools sort -o "${sample_id}.bam" ${mapped}
    """
}