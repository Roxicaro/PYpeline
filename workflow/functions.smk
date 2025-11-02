configfile: "config.yaml"

#Get input FASTQ files for FastQC pre-trimming
def get_fastqc_input_fastqs(wildcards):
    #Check if pair-end or single-read
    #Paired-end (Illumina)
    if config["seq_type"] == "PE":
        return [config["pe_samples"][wildcards.sample]["fastq1"],
                config["pe_samples"][wildcards.sample]["fastq2"]]
    #Single-read (Ion Torrent)
    elif config["seq_type"] == "SR":
        return config["fastq"]

#Get input FASTQ files for trimming
def get_trimming_input_fastqs(wildcards):
    #Check if pair-end or single-read
    #Paired-end (Illumina)
    if config["seq_type"] == "PE":
        return [config["pe_samples"][wildcards.sample]["fastq1"],
                config["pe_samples"][wildcards.sample]["fastq2"]]
    #Single-read (Ion Torrent)
    elif config["seq_type"] == "SR":
        return config["fastq"]

#Get input FASTQ files for BWA mapping
def get_bwa_map_input_fastqs(wildcards):
    #Check if pair-end or single-read
    #Paired-end (Illumina)
    if config["seq_type"] == "PE":
        if config["trim"] == "true":
            return [f"../results/trimmed_fastq/{wildcards.sample}_val_1.fq.gz", 
                    f"../results/trimmed_fastq/{wildcards.sample}_val_2.fq.gz"]  #get trimmed fastq files
        else:
            return [config["pe_samples"][wildcards.sample]["fastq1"],
                    config["pe_samples"][wildcards.sample]["fastq2"]] #get untrimmed fastq files
    #Single-read (Ion Torrent)
    elif config["seq_type"] == "SR":
        if config["trim"] == "true":
            return f"../results/trimmed_fastq/{wildcards.sample}_trimmed.fq.gz" #get trimmed fastq file
        else:
            return config["fastq"] #get untrimmed fastq file