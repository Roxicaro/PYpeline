nextflow.enable.dsl=2

//Processes block
process FASTQC {
    tag "$sample_id"

    conda "envs/fastqc.yaml"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}_fastqc")

    script:
    """
    mkdir ${sample_id}_fastqc
    fastqc -t ${task.cpus} -o ${sample_id}_fastqc ${reads.join(' ')}
    """
}

process TRIM_GALORE {
    tag "$sample_id"

    conda "envs/trim_galore.yaml"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz")

    script:
    if (params.tech == "illumina") {
        """
        trim_galore --fastqc --paired --gzip --output_dir . ${reads.join(' ')}
        """
    } else {
        """
        trim_galore --adapter ${params.adapter} --fastqc --gzip --output_dir . ${reads.join(' ')}
        """
    }
}

process BWA_MAP {
    tag "$sample_id"
    conda "envs/bwa.yaml"
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bwa mem -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' -t ${task.cpus} \
        ${params.ref} ${reads.join(' ')} \
        | samtools view -Sb - > ${sample_id}.bam
    """
}

process SAMTOOLS_SORT {
    tag "$sample_id"
    conda "envs/bwa.yaml"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools sort -O bam -o ${sample_id}.sorted.bam ${bam}
    """
}

process SAMTOOLS_INDEX {
    tag "$sample_id"
    conda "envs/bwa.yaml"

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sorted_bam}.bai")

    script:
    """
    samtools index ${sorted_bam}
    """
}

process MUTECT2 {
    tag "$sample_id"
    conda "envs/gatk.yaml"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path("${sample_id}.raw.vcf.gz")

    script:
    """
    gatk Mutect2 \
      -R ${params.ref} \
      -I ${bam} \
      -L ${params.bed} \
      -O ${sample_id}.raw.vcf.gz
    """
}

// Workflow block
workflow {
    /*
     * Dynamically build input channels based on seq_type (SR or PE).
     * Equivalent to your Snakemake logic using config["seq_type"].
     */
    Channel
        .fromFilePairs(params.reads, flat: true)
        .set { raw_reads }

    // PIPELINE (equivalent to Snakemake DAG)
    trimmed = TRIM_GALORE(raw_reads)
    qc      = FASTQC(raw_reads)
    mapped  = BWA_MAP(trimmed)
    sorted  = SAMTOOLS_SORT(mapped)
    indexed = SAMTOOLS_INDEX(sorted)
    called  = MUTECT2(indexed)

}