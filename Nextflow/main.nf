#!/usr/bin/env nextflow

// Include modules
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { BWA_ALIGN } from './modules/bwa_align.nf'
include { SAMTOOLS_SORTING } from './modules/samtools_sorting.nf'
include { SAMTOOLS_INDEXING } from './modules/samtools_indexing.nf'
include { GATK_MUTECT2 } from './modules/gatk_mutect2.nf'


workflow {
    raw_reads = Channel.fromFilePairs(params.reads)
    bwa_index = file(params.bwa_index)
    fasta = file(params.reference)
    dict = file(params.dict)

    FASTQC(raw_reads)                                       // Quality check of raw (untrimmed) FASTQ files
    trimmed_reads = TRIM_GALORE(raw_reads)                  // .fastq - Trimmed FASTQ files
    mapped = BWA_ALIGN(trimmed_reads, bwa_index)            // .bam - Unsorted aligned BAM files
    sorted = SAMTOOLS_SORTING(mapped)                       // .bam - Sorted BAM files
    indexed = SAMTOOLS_INDEXING(sorted)                     // .bai - BAI files (indexed BAM)
    GATK_MUTECT2 (sorted, indexed, bwa_index, fasta, dict)  // .raw.vcf.gz - Raw VCF files with variant calls
}