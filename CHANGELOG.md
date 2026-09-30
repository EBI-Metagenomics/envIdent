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

## v1.0.0 - 30.09.26

First version of the pipeline used to generate results for the EnvIdent database

### `Added`

Updated nextflow and nf-schema versions

Syntax changes to conform to the new schema

Updated the following:
 - DADA2
 - Cutadapt
 - Fastp
 - Seqfu
 - Seqtk
 - Seqkit
 - utils_nextflow_pipeflone subworkflow
 - MultiQC
 - dumpsoftwareversions
 - PIMENTO
 - Biopython

New option for users to supply their own primer sequences. New modules to reverse complement these primers and prepare them for cutadapt

New module to extract cutadapt prepared primers for the PIMENTO outputted primer sequences

Adapted the cutadapt module to run twice, with the original primers and then again with the reverse complemented primers

Added a mechanism to record the position at which primers were trimmed

Switched to running pimento for each strand

Updated the DADA2 subworkflow and DADA2.R, incorporating changes from our amplicon pipeline

Added the pimento bcv module

Created nf-tests for all modules missing one

Made an end to end test

Publishing the primers used by cutadapt

Support DADA2 merge-mode option

Publishing fastp jsons

Increased robustness of the pipeline to allow for low read count runs to be processed

Changed the reference databases to be optional inputs. Pipeline can be run without them or with one.

Added vsearch

Added LCA modules

Added vsearch formatting module

Added lca to krona count subworkflows and modules

Added tabix gzip module

Switch from using superkingdom to domain

Added download with fire module

Added run standardisation optional preprocessing step

Added BBMap module

Defined min percentage as a function - also allowing minimum to be set to 0

Added multiQC report per sample

Updated default parameters in modules.config

### `Fixed`

Bug fix - single-end reads were being run through seqkit seq twice. That's been resolved.

### `Deprecated`

Removed read merging prior to running pimento

Removed assessmcpproportions module

MapSeq has been replaced

Removal of mapseq to asv table module

Removed study level multiQC report
