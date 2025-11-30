#!/usr/bin/env nextflow

process TRIM_GALORE {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    conda "envs/trim_galore.yaml"
    container "community.wave.seqera.io/library/trim-galore_awscli:0051297f49a97842"
    //awscli        = 2.31.33   (conda-forge)
    //trim-galore   = 0.6.10    (bioconda)

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*_1_val_1.fq"), path("*_2_val_2.fq")

    script:
    """
    trim_galore --paired \
        --cores ${task.cpus} \
        --fastqc \
        --output_dir ./ \
        ${read1} ${read2}
    """
}