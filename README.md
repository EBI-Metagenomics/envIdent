

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

** Analysis:**
- Automatic  identification using PIMENTO
-  trimming using Cutadapt
-  validation and reporting

**Taxonomic Profiling:**
- Pfam-based COI (Cytochrome C Oxidase subunit I) profiling using HMMER
- Reads percentage threshold filtering for marker gene identification (configurable threshold)

**ASV Analysis:**
- Amplicon Sequence Variant (ASV) calling using DADA2
- ASV taxonomic classification using VSEARCH
- Krona chart visualization for taxonomic results

**Reporting and Quality Control:**
- Comprehensive MultiQC reports
- Failed and passed runs tracking
- Software version reporting

## Tools

| Tool | Version | Purpose |
|------|---------|---------|
| [cutadapt](https://cutadapt.readthedocs.io/en/stable/)  | 4.6 |  trimming |
| [DADA2](https://benjjneb.github.io/dada2/index.html)   | 1.30.0 | ASV calling and denoising |
| [fastp](https://github.com/OpenGene/fastp)  | 0.23.4 | Read quality control and filtering |
| [FastQC](https://github.com/s-andrews/fastqc) | 0.12.1 | Read quality control |
| [HMMER](http://hmmer.org/) | 3.4 | Profile HMM searching for COI sequences |
| [Krona](https://github.com/marbl/Krona)  | 2.8.1 | Interactive taxonomic visualization |
| [VSEARCH](https://github.com/torognes/vsearch)  | 2.32.0 | Taxonomic classification of ASVs |
| [mgnify-pipelines-toolkit](https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit) | 1.0.4 | Toolkit containing various in-house processing scripts |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 1.27 | Aggregated quality control reporting |
| [PIMENTO](https://github.com/EBI-Metagenomics/PIMENTO)  | 1.0.3 |  identification and inference |
| [SeqFu](https://telatin.github.io/seqfu2/) | 1.20.3 | FASTQ validity check |
| [SeqKit](https://bioinf.shenwei.me/seqkit/) | 2.9.0 | Read extraction and protein tranlsation |
| [Seqtk](https://github.com/lh3/seqtk) | 1.3 | Converting FASTQ to FASTA |

## Reference Databases

This pipeline uses the following reference databases:

| Database | Purpose | Default Location |
|----------|---------|------------------|
| BOLD | COI taxonomic classification and popular COI s | Configurable via parameters |
| MIDORI2 | COI taxonomic classification | Configurable via parameters |

> [!NOTE]
Running with both COI reference databases is enabled by default (`run_coi_bold = true` and `run_coi_midori = true`). Supply their `--coi_bold_ref_db` and `--coi_midori_ref_db` paths. Use `--run_coi_bold false` or `--run_coi_midori false` to skip a database. Set both run_coi_bold and run_coi_midori to false to generate ASVs and read counts without taxonomic assignments or Krona reports.
> Database paths can be configured in the pipeline parameters. Contact the development team for access to preprocessed databases.

## How to Run

### Requirements

The pipeline requires:
- Nextflow (≥24.04.2)
- Docker, Singularity, or Conda for software management
- Access to reference databases
- Database formatted for PIMENTO - a FASTA file with contig ids ending with F for forward strand and R for reverse strand. See [here](https://github.com/EBI-Metagenomics/PIMENTO/blob/main/pimento/standard_s/V3-V5.fasta) for an example

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
| `--std_primer_library` | `./data/standard_primers` | Directory containing forward (`*F.fasta`) and reverse (`*R.fasta`) PIMENTO primer libraries |
| `--cutadapt_primers` | `Default path` | Directory containing forward (`*F.fasta`), reverse (`*R.fasta`), reverse complemented forward (`F_RC.fasta`) and reverse complemented reverse (`R_RC.fasta`) cutadapt prepared primer libraries |
| `--pfam_coi_db` | `Default path` | Path to Pfam COI HMM database |


## Outputs

### Output directory structure

Example output structure for a sample (sample1). The qc_passed and qc_failed csvs are only present if you have samples that passed or failed:
```bash
results/
├── sample1/
│   ├── asv/
│   │   ├── sample1_asv_read_counts.tsv
│   │   ├── sample1_dada2_stats.tsv
│   │   └── sample1_asvs.fasta
│   ├── hmmsearch-COI/
│   │   ├── sample1_Pfam-A.domtbl
│   │   └── sample1_Pfam-A.txt
│   ├── primer-identification/
│   │   └── sample1.cutadapt.json
│   ├── qc/
│   │   ├── sample1_seqfu.tsv
│   │   └── sample1.fastp.json
│   │   └── sample1_1.fastq.gz
│   │   └── sample1_2.fastq.gz
│   │   └── sample1_suffix_header_err.json
│   ├── taxonomy-summary/
│   │   ├── BOLD/
│   │   │   ├── sample1_BOLD_vsearch_raw_hits.tsv
│   │   │   ├── sample1_BOLD_vsearch_hits_with_accessions.tsv
│   │   │   ├── sample1_BOLD_vsearch_hits_for_lca.tsv
│   │   │   ├── sample1_BOLD_taxonomy_lca_all_hits.tsv
│   │   │   ├── sample1_BOLD_taxonomy_lca_top_hits.tsv
│   │   │   ├── sample1_BOLD_krona_lca_all_hits_counts.tsv
│   │   │   ├── sample1_BOLD_krona_lca_top_hits_counts.tsv
│   │   │   ├── sample1_BOLD_krona_lca_all_hits.html
│   │   │   └── sample1_BOLD_krona_lca_top_hits.html
│   │   ├── MIDORI/
│   │   │   ├── sample1_MIDORI_vsearch_raw_hits.tsv
│   │   │   ├── sample1_MIDORI_vsearch_hits_with_accessions.tsv
│   │   │   ├── sample1_MIDORI_vsearch_hits_for_lca.tsv
│   │   │   ├── sample1_MIDORI_taxonomy_lca_all_hits.tsv
│   │   │   ├── sample1_MIDORI_taxonomy_lca_top_hits.tsv
│   │   │   ├── sample1_MIDORI_krona_lca_all_hits_counts.tsv
│   │   │   ├── sample1_MIDORI_krona_lca_top_hits_counts.tsv
│   │   │   ├── sample1_MIDORI_krona_lca_all_hits.html
│   │   │   └── sample1_MIDORI_krona_lca_top_hits.html
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

Taxonomy filenames use `<sample>_<database>_<description>`. The database labels
come from `--bold_label` (default `BOLD`) and `--midori_label`
(default `MIDORI`); each label controls both the folder and filename prefix.

| Description | Contents |
| --- | --- |
| `vsearch_raw_hits.tsv` | Original VSEARCH output, without a header |
| `vsearch_hits_with_accessions.tsv` | Readable hits with the standard taxonomy ranks, retaining reference database accessions |
| `vsearch_hits_for_lca.tsv` | Individual hits with taxonomy, identity and query coverage for LCA script input |
| `taxonomy_lca_all_hits.tsv` | ASV assignments calculated from all hits above the percentage identity threshold at the selected rank |
| `taxonomy_lca_top_hits.tsv` | ASV assignments from highest-identity qualifying hits, above the percentage identity threshold at the selected rank |
| `krona_lca_all_hits_counts.tsv`, `krona_lca_top_hits_counts.tsv` | Headerless, read-weighted counts grouped by taxonomy; includes unclassified reads |
| `krona_lca_all_hits.html`, `krona_lca_top_hits.html` | Interactive reports for the respective assignment method |

For LCA input, the matching genus prefix is removed from the species label first.
If the first remaining underscore-separated word contains `.`, eg. `sp.`, the species rank
is left empty (`s__;`). Hits in the accessions output file retain the full label after genus removal.

Taxonomy uses eight ranks: domain, kingdom, phylum, class, order,
family, genus, and species. Missing ranks retain empty placeholders.

ASVs without reference hits remain in the output tables. Raw VSEARCH rows use
`*` for the missing target. Both formatted tables and LCA assignment tables
retain the ASV ID with `d__;k__;p__;c__;o__;f__;g__;s__;`. Formatted no-hit rows
have an empty accession and `NA` identity/coverage where those columns are present.
Their reads remain included under `Unclassified` in the Krona counts and reports.

`asv/<sample>_asv_read_counts.tsv` is generated once per sample by
the MGnify pipelines toolkit's `make_asv_count_table`, counting nonzero forward-map
entries for the filtered reads from DADA2. Counts are independent of taxonomy and
are generated even when both database branches are disabled. BOLD and MIDORI
share this table. Empty counts retain the `asv` and `count` header.

Both final `taxonomy_lca_all_hits.tsv` and `taxonomy_lca_top_hits.tsv` files are
headerless, with columns `ASV ID`, `taxonomy`, and `count`. The third column comes
from the shared ASV read counts table.

### Key Output Files

* **MultiQC Report**: Comprehensive quality control summary across all samples
* **ASV Sequences**: FASTA files containing called Amplicon Sequence Variants
* **Taxonomic Classifications**: TSV files with taxonomic assignments for ASVs
* **Krona Charts**: Interactive HTML visualisations of taxonomic composition
* **QC Summary Files**: Lists of samples that passed or failed quality control steps

### Configuration Profiles

`--dada2_merge_mode` supports `standard` and `gap`. Separate-strand mode is
currently disabled because ASV counting does not support its `seq_f_N` and
`seq_r_N` identifiers.

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
