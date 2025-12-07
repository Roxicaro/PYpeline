# PYpeline
**A lightweight and modular NGS pipeline from FASTQ → BAM → VCF**, designed for small projects and rapid iteration.\
Ready for AWS Batch.

---

## Key Features

### Implemented in two workflow languages
- **Snakemake version** (ideal for development, rule-based execution locally)
- **Nextflow version** (designed for cloud execution and scalability)

---

This pipeline automates the following steps:
- Trimming (TrimGalore)
- Quality control (FastQC + MultiQC)
- Read alignment (BWA MEM)
- Sorting and indexing (Samtools)
- Variant calling (Mutect2)

It supports **Single-End (IonTorrent, AmpliSeq)** (_Snakemake only_) and **Paired-End (Illumina)** FASTQ files.

## Quickstart (TL;DR)

- Clone the repo:
```bash
git clone https://github.com/Roxicaro/Pypeline.git
cd Pypeline
```
- Install requirements: Snakemake / Nextlow (Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through WSL.)
- Prepare `data/` directory (FASTQs + reference)
- Edit `workflow/config.yaml` (Snakemake) or `nexflow.config` (Nexflow)
- Add Sample_ID and fastq file paths to data/samplesheet.csv (Nextflow)
- Run


## Directory Structure
Before running the pipeline, organize your data as follows:

### Snakemake
```markdown
data/
├── fastq_files/      # FASTQ input files
├── bed/              # BED files for target regions (optional)
└── references/       # Reference genome files (FASTA + indexe files for BWA MEM and Mutect2)

workflow/
├── envs/             # Environment .yaml files
├── config.yaml       # Configuration file where the user sets pipeline parameters and file paths
├── Snakefile         # Main workflow file              
└── functions.smk     # Python functions to get input files

results/              # Output files will be written here
```

### Nextflow
```markdown
Nextflow/
├── envs/                 # Environment .yaml files
├── main.nf               # Main workflow file              
├── nextflow.config       # Configuration file where the user sets pipeline parameters and file paths
├── results/              # Output files will be written here
├── modules/              # Code for the different workflow modules
└── data/
    ├── fastq_files/      # FASTQ input files
    ├── samplesheet.csv   # Sample list and file paths
    ├── bed/              # BED files for target regions (optional)
    └── references/       # Reference genome files (FASTA + indexe files for BWA MEM and Mutect2)
```

- **FASTQ files:** Input sequencing reads.  
- **BED files:** Target regions for variant calling (optional).  
- **References:** Reference genome files including any required index files for `bwa mem` and `Mutect2` _(.fa / .fai / .dict / .amb / .ann / .bwt / .pac / .sa)_.

## Input Sample Sheet (Required for Nextflow)

To run the Nexflow pipeline, you must provide a **samplesheet CSV** listing the input FASTQ files.\
This file **must be located at**: `Nextflow/data/samplesheet.csv`

### CSV Format

| sample_id | read1 | read2 |
|----------|------|-------|
| Unique sample name | Path to R1 FASTQ | Path to R2 FASTQ |

---

### Example Sample Sheet

```csv
sample_id,read1,read2
SRR35855706,data/fastq_files/SRR35855706_1.fastq,data/fastq_files/SRR35855706_2.fastq
```

## Running the pipeline
First, edit `config.yaml` (Snakemake) or `data/samplesheet.csv` (Nextflow) to **inform FASTQ file paths**.\
To run using AWS Batch (Nextflow) edit `nextflow.conf` to include the paths to the reference and index files in an S3 bucket.

### Running with Snakemake (local w/ Docker):
```markdown
docker run -it --rm \
  -v $PWD:/pipeline \
  -w /pipeline/workflow \
  roxicaro/pypeline-snakemake \
  snakemake --cores 4
```
`--cores` specifies the number of cores to be used.

### Running with Nextflow:
**Local (Not recommended. Requires all tools to be locally installed):**
```markdown
nextflow run main.nf
```
**Docker (Recommended):**
```markdown
nextflow run main.nf -profile docker
```
**AWS Batch:**
```markdown
nextflow run main.nf -profile aws_batch
```

### Generate a DAG diagram showing the workflow (Snakemake):
```markdown
docker run -it --rm \
  -v $PWD:/pipeline \
  -w /pipeline/workflow \
  roxicaro/pypeline-snakemake \
  snakemake -np --dag | dot -Tsvg > dag.svg
```

## Output structure
After processing, results are written to:
```markdown
results/
├── trimmed_fastq/    # FASTQ files after trimming
├── fastqc/           # QC reports (per sample + MultiQC)
├── mapped_reads/     # Unsorted BAM
├── sorted_reads/     # Sorted BAM and BAI files
├── variant_calls/    # `.vcf` files
└── logs/             # Execution logs for troubleshooting
```

## Requirements
- [Snakemake](https://snakemake.github.io/) / [Nextflow](https://github.com/nextflow-io/nextflow) \
Note: Nextflow can be used on any POSIX-compatible system (Linux, macOS, etc), and on Windows through WSL.

## Example Workflow DAG (Snakemake)

Pair-end sequencing FASTQs (Illumina):

![Pipeline DAG](Snakemake/workflow/dag_pe.svg)

Single-read sequencing FASTQs (Ion Torrent):

![Pipeline DAG](Snakemake/workflow/dag_sr.svg)

## Example Workflow DAG (Nextflow)
![Pipeline DAG](Nextflow/dag.svg)
