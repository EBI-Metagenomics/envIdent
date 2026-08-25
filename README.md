

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# EnvIdent - EBI-Metagenomics eDNA Analysis Pipeline

This repository contains EnvIdent v0.1 - EBI-Metagenomic's eDNA analysis pipeline. This pipeline is designed for the analysis of environmental DNA (eDNA) sequencing data, implementing a comprehensive workflow for quality control, primer identification, Amplicon Sequence Variant (ASV) calling and taxonomic profiling using modern bioinformatics tools.

Currently the pipeline supports analysis of COI metabarcoding reads.

## Pipeline Description

<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/envident_schema.png">
    <img alt="EBI-Metagenomics/envident" src="docs/images/envident_schema.png" style="width: 80%;">
  </picture>
</h1>


### Features

EnvIdent v0.1 implements the following key features:

**Quality Control and Preprocessing:**
- Raw reads quality assessment using FastQC
- Reads quality control and filtering using fastp
- Minimum read count filtering (configurable threshold)

**Primer Analysis:**
- Automatic primer identification using PIMENTO
- Primer trimming using Cutadapt
- Primer validation and reporting

**Taxonomic Profiling:**
- Pfam-based COI (Cytochrome C Oxidase subunit I) profiling using HMMER
- Reads percentage threshold filtering for marker gene identification (configurable threshold)

**ASV Analysis:**
- Amplicon Sequence Variant (ASV) calling using DADA2
- ASV taxonomic classification using MAPseq
- Krona chart visualization for taxonomic results

**Reporting and Quality Control:**
- Comprehensive MultiQC reports
- Failed and passed runs tracking
- Software version reporting

## Tools

| Tool | Version | Purpose |
|------|---------|---------|
| [cutadapt](https://cutadapt.readthedocs.io/en/stable/)  | 4.6 | Primer trimming |
| [DADA2](https://benjjneb.github.io/dada2/index.html)   | 1.30.0 | ASV calling and denoising |
| [fastp](https://github.com/OpenGene/fastp)  | 0.23.4 | Read quality control and filtering |
| [FastQC](https://github.com/s-andrews/fastqc) | 0.12.1 | Read quality control |
| [HMMER](http://hmmer.org/) | 3.4 | Profile HMM searching for COI sequences |
| [Krona](https://github.com/marbl/Krona)  | 2.8.1 | Interactive taxonomic visualization |
| [MAPseq](https://github.com/meringlab/MAPseq)  | 2.1.1b | Taxonomic classification of ASVs |
| [mgnify-pipelines-toolkit](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit) | 1.0.4 | Toolkit containing various in-house processing scripts |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 1.27 | Aggregated quality control reporting |
| [PIMENTO](https://github.com/EBI-Metagenomics/PIMENTO)  | 1.0.3 | Primer identification and inference |
| [SeqFu](https://telatin.github.io/seqfu2/) | 1.20.3 | FASTQ validity check |
| [SeqKit](https://bioinf.shenwei.me/seqkit/) | 2.9.0 | Read extraction and protein tranlsation |
| [Seqtk](https://github.com/lh3/seqtk) | 1.3 | Converting FASTQ to FASTA |

## Reference Databases

This pipeline uses the following reference databases:

| Database | Purpose | Default Location |
|----------|---------|------------------|
| BOLD | COI taxonomic classification and popular COI primers | Configurable via parameters |
| MIDORI2 | COI taxonomic classification | Configurable via parameters |

> [!NOTE]
> Database paths can be configured in the pipeline parameters. Contact the development team for access to preprocessed databases.

## How to Run

### Requirements

The pipeline requires:
- Nextflow (≥24.04.2)
- Docker, Singularity, or Conda for software management
- Access to reference databases
- Primer database formatted for PIMENTO - a FASTA file with contig ids ending with F for forward primers and R for reverse primers. See [here](https://github.com/EBI-Metagenomics/PIMENTO/blob/main/pimento/standard_primers/V3-V5.fasta) for an example

### Input Format

The input data should be eDNA sequencing reads (paired-end or single-end) in FASTQ format, specified using a CSV samplesheet:

```csv
sample,fastq_1,fastq_2,single_end
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,false
sample2,/path/to/sample2.fastq.gz,,true
```
> [!NOTE]
> EnvIdent has not yet been optimised for single-end reads, the parameters used as default may not be optimal

### Basic execution

```bash
nextflow run EBI-Metagenomics/envident \
    -r main \
    -profile example_slurm \
    --input samplesheet.csv \
    --outdir results
```

### Key Parameters

| **Parameter** | **Default** | **Description** |
|----------------|-------------|-----------------|
| `--min_read_count` | `5000` | Minimum number of reads required per sample |
| `--reads_percentage_threshold` | `0.10` | Minimum percentage of reads matching COI profile |
| `--std_primer_library_forward` | `[]` | Forward primer library for PIMENTO R1 identification |
| `--std_primer_library_reverse` | `[]` | Reverse primer library for PIMENTO R2 identification |
| `--pfam_coi_db` | `Default path` | Path to Pfam COI HMM database |


## Outputs

### Output directory structure

Example output structure for a sample (sample1). The qc_passed and qc_failed csvs are only present if you have samples that passed or failed:
```bash
results/
├── sample1/
│   ├── asv/
│   │   ├── sample1_DADA2-BOLD_asv_read_counts.tsv
│   │   ├── sample1_DADA2-MIDORI_asv_read_counts.tsv
│   │   └── sample1_dada2_stats.tsv
│   │   └── sample1_asvs.fasta
│   ├── hmmsearch-COI/
│   │   ├── sample1_Pfam-A.domtbl
│   │   └── sample1_Pfam-A.txt
│   ├── primer-identification/
│   │   └── sample1.cutadapt.json
│   ├── qc/
│   │   ├── sample1_seqfu.tsv
│   │   └── sample1.fastp.json
│   │   └── sample1.merged.fastq.gz
│   │   └── sample1_suffix_header_err.json
│   ├── taxonomy-summary/
│   │   ├── DADA2-BOLD/
│   │   |   ├── ERR8441464_DADA2-BOLD_asv_krona_counts.txt
│   │   |   ├── ERR8441464_DADA2-BOLD_asv_taxa.tsv
│   │   |   ├── ERR8441464_DADA2-BOLD.html
│   │   |   └── ERR8441464_DADA2-BOLD.mseq
│   │   ├── DADA2-MIDORI/
│   │   |   ├── ERR8441464_DADA2-MIDORI_asv_krona_counts.txt
│   │   |   ├── ERR8441464_DADA2-MIDORI_asv_taxa.tsv
│   │   |   ├── ERR8441464_DADA2-MIDORI.html
│   │   |   └── ERR8441464_DADA2-MIDORI.mseq
├── pipeline_info/
│   ├── execution_report_YYYY-MM-DD_HH-mm-ss.html
│   ├── execution_timeline_YYYY-MM-DD_HH-mm-ss.html
│   ├── execution_trace_YYYY-MM-DD_HH-mm-ss.txt
│   ├── params_YYYY-MM-DD_HH-mm-ss.json
│   ├── pipeline_dag_YYYY-MM-DD_HH-mm-ss.html
│   └── envident_software_mqc_versions.yml
├── multiqc_report.html
├── qc_passed_runs.csv
└── qc_failed_runs.csv
```

### Key Output Files

* **MultiQC Report**: Comprehensive quality control summary across all samples
* **ASV Sequences**: FASTA files containing called Amplicon Sequence Variants
* **Taxonomic Classifications**: TSV files with taxonomic assignments for ASVs
* **Krona Charts**: Interactive HTML visualisations of taxonomic composition
* **QC Summary Files**: Lists of samples that passed or failed quality control steps

### Configuration Profiles

The pipeline includes pre-configured profiles:

* docker: Use Docker containers
* singularity: Use Singularity containers
* conda: Use Conda environments
* example_slurm: Optimized for SLURM clusters
* test: Small test dataset for validation

## Citations

This pipeline uses code developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
