configfile: "config.yaml"


#Aligning FASTQ reads with BWA MEM and generate SAM/BAM output with samtools
def get_bwa_map_input_fastqs(wildcards):
    #Check if pair-end or single-read
    #Paired-end (Illumina)
    if config["seq_type"] == "PE":
        return config["fastq1"], config["fastq2"]
    #Single-read (Ion Torrent)
    elif config["seq_type"] == "SR":
        return config["fastq"]