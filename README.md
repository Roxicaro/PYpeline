# PYpeline
Basic NGS pipeline from FASTQ to vcf

## Directory Structure
Before running the pipeline, organize your data as follows:

```markdown
data/
├── fastq_files/      # FASTQ input files
├── bed/              # BED files for target regions
└── references/       # Reference genome files (FASTA + indexes for BWA MEM and Mutect2)
```

- **FASTQ files:** Input sequencing reads.  
- **BED files:** Target regions for variant calling.  
- **References:** Reference genome files including any required index files for `bwa mem` and `Mutect2`.

## Workflow DAG

Pair-end sequencing FASTQs (Illumina):
![Pipeline DAG](workflow/dag_pe.svg)

Single-read sequencing FASTQs (Ion Torrent):
![Pipeline DAG](workflow/dag_pe.svg)