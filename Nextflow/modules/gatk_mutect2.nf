#!/usr/bin/env nextflow

process GATK_MUTECT2 {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calls", mode: 'copy'

    conda "envs/gatk.yaml"
    container "community.wave.seqera.io/library/gatk4_awscli:ecc06537e2369f27"
    //awscli	= 2.31.33   (conda-forge)
    //gatk4     = 4.6.2.0   (bioconda)

    input:
    tuple val(sample_id), path(sorted)
    path(indexed)
    path bwa_index
    path(fasta)
    path(dict)

    output:
    path "${sample_id}.raw.vcf.gz"

    script:
    """
    gatk Mutect2 -R ${fasta} \
    -I ${sorted} \
    -O "${sample_id}.raw.vcf.gz"
    """
}