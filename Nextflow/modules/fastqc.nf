#!/usr/bin/env nextflow

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    conda "envs/fastqc.yaml"
    container "community.wave.seqera.io/library/fastqc_awscli:6c9886bd361f561a"
    //awscli    = 2.31.33   (conda-forge)
    //fastqc    = 0.12.1    (bioconda)

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2}
    """
}