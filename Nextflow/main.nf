#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    conda "envs/fastqc.yaml"

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
        --output_dir ./ \
        ${reads.join(" ")}
    """
}

workflow {
    Channel.fromFilePairs(params.reads)
        .set { raw_reads }

    FASTQC(raw_reads)
    TRIM_GALORE(raw_reads)
}