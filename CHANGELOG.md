# EBI-Metagenomics/envident: Changelog

## v0.0.1 - 21.05.26

Initial release of EBI-Metagenomics/envident.

### `Added`

Added support for hmmsearch to use read counts
Added params to fastp
Updated pimento to latest version, using multiple threads
Adjusted minimum reads percentage threshold
Readme, configs and json files updated
Updated hmmsearch and seqkit modules
Renamed to EnvIdent

### `Fixed`

Bug fix for primer splitting
Bug fix to ensure cutadapt is run on the correct reads

### `Deprecated`

Removed igenomes
Removed unused modules and subworkflows


## v0.0.2 - 19.08.26

Syntax updates required for nextflow, nf-schema and nf-prov updates

### `Fixed`

Using _ isn't allowed as an identifier anymore, renamed
Nextflow strict syntax doesn't support import declarations
Errors due to changes in how variables are defined
Re-use of variable name in the surrounding scope
Issues with multiqc file publishing
Preventing warning with ref db definition
Dada2 container path to allow working with docker
