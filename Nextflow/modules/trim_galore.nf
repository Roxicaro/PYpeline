#!/usr/bin/env nextflow

process TRIM_GALORE {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    conda "envs/trim_galore.yaml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_val_1.fq"), path("${sample_id}_2_val_2.fq")

    script:
    """
    trim_galore --paired \
        --cores ${task.cpus} \
        --fastqc \
        --output_dir ./ \
        ${reads.join(" ")}
    """
}