#!/usr/bin/env nextflow

process BWA_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned_unsorted", mode: 'copy'

    conda "envs/bwa.yaml"
    container "community.wave.seqera.io/library/bwa_samtools_awscli:8ec39ef35cae626f"
    //awscli	= 2.31.33   (conda-forge)
    //bwa	    = 0.7.19    (bioconda)
    //samtools  = 1.22.1    (bioconda)

    input:
    tuple val(sample_id), path(read1), path(read2)
    path fasta
    path bwa_index  //stages reference indexes

    output:
    tuple val(sample_id), path("${sample_id}_unsorted.bam")

    script:
    def param = "@RG\\tID:${sample_id}\\tSM:${sample_id}" // read group parameter
    """
    bwa mem -R '${param}' -t ${task.cpus} ${fasta} ${read1} ${read2} | samtools view -b - > "${sample_id}_unsorted.bam"
    """
}