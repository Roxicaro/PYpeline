configfile: "config.yaml"

#Get input FASTQ files for FastQC pre-trimming
def get_preprocessed_input_fastqs(wildcards):
    #Check if pair-end or single-read
    #Paired-end (Illumina)
    if config["seq_type"] == "PE":
        return [config["pe_samples"][wildcards.sample]["fastq1"],
                config["pe_samples"][wildcards.sample]["fastq2"]]
    #Single-read (Ion Torrent)
    elif config["seq_type"] == "SR":
        return config["sr_samples"][wildcards.sample]

#Get input FASTQ files for BWA mapping
def get_bwa_map_input_fastqs(wildcards):
    if config["trim"] == "true":
        if config["seq_type"] == "PE":
            return [f"../results/trimmed_fastq/{wildcards.sample}_1_val_1.fq.gz",
                    f"../results/trimmed_fastq/{wildcards.sample}_2_val_2.fq.gz"]
        else:
            return f"../results/trimmed_fastq/{wildcards.sample}_trimmed.fq.gz" 
    else:
        return get_trimming_input_fastqs(wildcards)