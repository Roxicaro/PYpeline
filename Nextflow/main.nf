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
        --fastqc \
        --output_dir ./ \
        ${reads.join(" ")}
    """
}

process BWA_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned_unsorted", mode: 'copy'

    conda "envs/bwa.yaml"

    input:
    tuple val(sample_id), path(read1), path(read2)
    path bwa_index

    output:
    tuple val(sample_id), path("${sample_id}_unsorted.bam")

    script:
    def idxbase = bwa_index[0].baseName // sets the index base name
    def param = "@RG\\tID:${sample_id}\\tSM:${sample_id}" // read group parameter
    """
    bwa mem -R '${param}' -t ${task.cpus} ${idxbase} ${read1} ${read2} | samtools view -b - > "${sample_id}_unsorted.bam"
    """
}

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

process GATK_MUTECT2 {
    tag "$sample_id"
    publishDir "${params.outdir}/variant_calls", mode: 'copy'

    conda "envs/gatk.yaml"

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
workflow {
    raw_reads = Channel.fromFilePairs(params.reads)
    //bwa_index = file( 'data/references/hg38.fa.{,amb,ann,bwt,pac,sa,fai}' )
    bwa_index = file(params.bwa_index)
    fasta = file(params.reference)
    dict = file(params.dict)

    FASTQC(raw_reads)
    trimmed_reads = TRIM_GALORE(raw_reads)
    mapped = BWA_ALIGN(trimmed_reads, bwa_index)
    sorted = SAMTOOLS_SORTING(mapped)
    indexed = SAMTOOLS_INDEXING(sorted)
    GATK_MUTECT2 (sorted, indexed, bwa_index, fasta, dict)
}