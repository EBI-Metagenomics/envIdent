/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_QC                     } from '../subworkflows/local/reads_qc/main'
include { READS_QC as READS_QC_BEFOREHMM   } from '../subworkflows/local/reads_qc/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FROM_FIRE } from '../modules/ebi-metagenomics/downloadfromfire/main'
include { BBMAP_REFORMAT_STANDARDISE } from '../modules/ebi-metagenomics/bbmap/reformat_standardise/main'
include { FASTQC as FASTQC_RAW         } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CLEAN       } from '../modules/nf-core/fastqc/main'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_nfcore_envident_pipeline'
include { PRIMER_IDENTIFICATION as PRIMER_IDENTIFICATION_F } from '../subworkflows/local/primer_identification_swf.nf'
include { PRIMER_IDENTIFICATION as PRIMER_IDENTIFICATION_R } from '../subworkflows/local/primer_identification_swf.nf'
include { CONCAT_PRIMER_CUTADAPT       } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { SUPPLIED_PRIMERS             } from '../modules/local/supplied_primers/main.nf'
include { PREP_CUTADAPT_PRIMERS        } from '../modules/local/prep_cutadapt_primers/main.nf'
include { EXTRACT_CUTADAPT_PRIMERS     } from '../modules/local/extract_cutadapt_primers/main.nf'
include { PROFILE_HMMSEARCH_PFAM       } from '../subworkflows/local/profile_hmmsearch_pfam/main'
include { DADA2_SWF                    } from '../subworkflows/local/dada2_swf.nf'
include { VSEARCH_ASV_LCA as VSEARCH_ASV_LCA_BOLD } from '../subworkflows/local/vsearch_asv_krona_swf/main.nf'
include { VSEARCH_ASV_LCA as VSEARCH_ASV_LCA_MIDORI2 } from '../subworkflows/local/vsearch_asv_krona_swf/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../modules/local/make_asv_count_tables/main'
include { LCA_KRONA_REPORTS as LCA_KRONA_REPORTS_BOLD } from '../subworkflows/local/vsearch_asv_krona_swf/main'
include { LCA_KRONA_REPORTS as LCA_KRONA_REPORTS_MIDORI2 } from '../subworkflows/local/vsearch_asv_krona_swf/main'
include { MULTIQC                      } from '../modules/nf-core/multiqc/main'

// Import samplesheetToList from nf-schema //
include { samplesheetToList            } from 'plugin/nf-schema'

// Import reads_merged_input_prep function (it's very big and deserved to be in its own file) //
include { reads_merged_input_prep      } from '../bin/reads_merged_input_prep.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ENVIDENT {

    take:
    samplesheet // channel: samplesheet read in from --input
    main:
    
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INITIALISE PRIMER LIBRARIES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Initialiase standard primer library for PIMENTO if user-given//
    // If there are no primers provided, it will fallback to use the default PIMENTO standard primer library
    std_primer_library = []

    if (params.std_primer_library){
        std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)
    }

    cutadapt_primers = []

    if (params.cutadapt_primers){
        cutadapt_primers = file(params.cutadapt_primers, type: 'dir', checkIfExists: true)
    }

    // Organise input tuple channel //
    def groupReads = { meta, fq1, fq2 ->
        def single_file = (fq2 == [])
        meta['interleaved'] = (!meta.single_end) && single_file
        if (single_file) {
            return tuple(meta, [fq1])
        }
        else {
            return tuple(meta, [fq1, fq2])
        }
    }
    ch_input = samplesheet.map(groupReads)

    if (params.use_fire_download) {
        /*
         * Internally we need to bypass Nextflow S3 integration until https://github.com/nextflow-io/nextflow/issues/4873 is fixed
         * The EBI parameter is needed as this only works on EBI network, FIRE is not accessible otherwise
        */
        DOWNLOAD_FROM_FIRE(
            ch_input
        )

        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE.out.versions.first())
        ch_input = DOWNLOAD_FROM_FIRE.out.downloaded_files
    }
    
    // Standardise headers and de-interleave as needed
    if (!params.skip_standardise) {
        standardise_input = ch_input.multiMap{
            meta, reads ->
            reads: [meta, reads]
            interleaved: meta.interleaved
        }
        BBMAP_REFORMAT_STANDARDISE(
            standardise_input.reads, 
            standardise_input.interleaved, 
            'fastq.gz'
        )
        ch_versions = ch_versions.mix(BBMAP_REFORMAT_STANDARDISE.out.versions)
        ch_input = BBMAP_REFORMAT_STANDARDISE.out.reformated
    }
    
    FASTQC_RAW(
        ch_input.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.id = meta.id + "_raw"
            [new_meta, reads]
        }
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    // Sanity checking and quality control of reads //
    READS_QC(
        true, 
        ch_input,
        false,
        false
    )
    ch_versions = ch_versions.mix(READS_QC.out.versions)

    // Filter and branch reads based on minimum read count with logging
    READS_QC.out.reads.branch{ meta, reads ->
                                    def read_files = reads instanceof List ? reads : [reads]
                                    def count = read_files[0].toAbsolutePath().countFastq()
                                    qc_pass: count >= params.min_read_count
                                    qc_fail: count < params.min_read_count
                                }
                                .set { extended_reads_qc }

    supplied_primers = extended_reads_qc.qc_pass
        .filter { meta, _reads -> meta.forward_primer && meta.reverse_primer }
        .map { meta, _reads ->
            meta + [direction: 'provided', direction_size: 0]
        }

    // Generates fasta files with the user supplied primer sequences from the samplesheet. Skipped if no primers provided
    SUPPLIED_PRIMERS(supplied_primers)

    // Takes user supplied primers and prepares them for cutadapt, including reverse complementing and adding the syntax needed for internal adaptors //
    PREP_CUTADAPT_PRIMERS(
        SUPPLIED_PRIMERS.out.supplied_primer_out
    )
    ch_versions = ch_versions.mix(PREP_CUTADAPT_PRIMERS.out.versions.first())

    primers_to_identify = extended_reads_qc.qc_pass
        .filter { meta, _reads -> !(meta.forward_primer && meta.reverse_primer) }

    primer_reads_f = primers_to_identify
        .map { meta, reads ->
            tuple(meta + [direction: 'f'], reads instanceof List ? reads[0] : reads)
        }

    primer_reads_r = primers_to_identify
        .map { meta, reads ->
            tuple(meta + [direction: 'r'], reads instanceof List ? reads[1] : reads)
    }

    // Identify primers independently in read 1 and read 2 for samples without supplied primers.
    PRIMER_IDENTIFICATION_F(
        primer_reads_f,
        std_primer_library
    )
    ch_versions = ch_versions.mix(PRIMER_IDENTIFICATION_F.out.versions)

    PRIMER_IDENTIFICATION_R(
        primer_reads_r,
        std_primer_library
    )
    ch_versions = ch_versions.mix(PRIMER_IDENTIFICATION_R.out.versions)   

    primer_outputs = PRIMER_IDENTIFICATION_F.out.std_primer_out
    .map { meta, primers -> [meta.subMap('id', 'single_end'), meta, primers] }
    .join(
        PRIMER_IDENTIFICATION_R.out.std_primer_out
            .map { meta, primers -> [meta.subMap('id', 'single_end'), primers] },
        by: [0]
    )
    .map { _key, meta, f_primers, r_primers ->
        tuple(meta + [direction: 'identified', direction_size: 0], f_primers, r_primers)
    }

    // Using the contig id of the primers identified by PIMENTO, extracts the corresponding formatted and revcomped primers from the cutadapt primer library for use in cutadapt trimming
    EXTRACT_CUTADAPT_PRIMERS(
        primer_outputs,
        cutadapt_primers
    )
    ch_versions = ch_versions.mix(EXTRACT_CUTADAPT_PRIMERS.out.versions.first())

    // Concatenate all primers for for a run, using either the user supplied primers or those identified by PIMENTO, send them to cutadapt with original QCd reads for primer trimming
    CONCAT_PRIMER_CUTADAPT(
        EXTRACT_CUTADAPT_PRIMERS.out.cutadapt_prepared_primers
            .mix(PREP_CUTADAPT_PRIMERS.out.cutadapt_prepared_primers),
        READS_QC.out.reads
    )
    ch_versions = ch_versions.mix(CONCAT_PRIMER_CUTADAPT.out.versions)

    reads_merge_input = reads_merged_input_prep(READS_QC.out.reads, CONCAT_PRIMER_CUTADAPT.out.cutadapt_out)

    FASTQC_CLEAN(
        reads_merge_input.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.id = meta.id + "_clean"
            [new_meta, reads]
        }
    )
    ch_versions = ch_versions.mix(FASTQC_CLEAN.out.versions.first())

    READS_QC_BEFOREHMM(
        false, 
        reads_merge_input,
        false,
        true // merge
    )
    ch_versions = ch_versions.mix(READS_QC_BEFOREHMM.out.versions) 

    // Exclude samples with no FASTA records after QC/merging from HMM profiling.
    reads_before_hmm = READS_QC_BEFOREHMM.out.reads_fasta
        .branch { _meta, reads ->
            qc_pass: reads.countFasta() > 0
            qc_empty: true
        }

    // Pfam profiling
    pfam_db = params.pfam_coi_db ?
    channel
        .fromPath(params.pfam_coi_db, checkIfExists: true)
        .first() :
    channel.empty()

    PROFILE_HMMSEARCH_PFAM(
        reads_before_hmm.qc_pass,
        pfam_db,
        READS_QC_BEFOREHMM.out.fastp_summary_json
    )
    ch_versions = ch_versions.mix(PROFILE_HMMSEARCH_PFAM.out.versions)
        
    // Filter samples based on reads_percentage threshold and get filtered domtbl
    ch_passed_samples = PROFILE_HMMSEARCH_PFAM.out.profile
        .filter { meta, tsv_file ->
            def threshold = params.reads_percentage_threshold ?: 0.10
            
            try {
                def lines = tsv_file.readLines()
                def dataLine = lines[1] // Skip header, get the single data row
                def columns = dataLine.split('\t')
                def readsPercentageStr = columns[4]
                
                def readsPercentage = readsPercentageStr as Double
                def passes = readsPercentage >= threshold
                
                return passes
                
            } catch (Exception e) {
                return false
            }
        }
        .map { meta, tsv_file -> meta }

    ch_filtered_domtbl = ch_passed_samples.join(PROFILE_HMMSEARCH_PFAM.out.domtbl) 

    // Run DADA2 ASV generation with filtered samples //
    DADA2_SWF(
        reads_merge_input,
        ch_filtered_domtbl
    )
    ch_versions = ch_versions.mix(DADA2_SWF.out.versions)
 
    def dada2_stats_fail = DADA2_SWF.out.dada2_stats_fail.map { meta, stats_fail ->
                                def key = meta.subMap('id', 'single_end')
                                return [key, stats_fail]
                            }

    // ASV taxonomic assignments + generate Krona plots for each run+amp_region //

    if (params.COI_bold_ref_db) {
        ref_db = file(params.COI_bold_ref_db, type: 'file', checkIfExists: true)
        VSEARCH_ASV_LCA_BOLD(
            DADA2_SWF.out.dada2_out,
            ref_db
        )
        ch_versions = ch_versions.mix(VSEARCH_ASV_LCA_BOLD.out.versions)
    }

    if (params.COI_midori2_ref_db) {
        ref_db = file(params.COI_midori2_ref_db, type: 'file', checkIfExists: true)
        VSEARCH_ASV_LCA_MIDORI2(
            DADA2_SWF.out.dada2_out,
            ref_db
        )
        ch_versions = ch_versions.mix(VSEARCH_ASV_LCA_MIDORI2.out.versions)
    }

    if (params.COI_bold_ref_db || params.COI_midori2_ref_db) {
        // Either complete LCA table supplies the ASV IDs, including no-hit ASVs.
        count_taxonomy = params.COI_bold_ref_db
            ? VSEARCH_ASV_LCA_BOLD.out.lca_all
            : VSEARCH_ASV_LCA_MIDORI2.out.lca_all
        map_count_input = DADA2_SWF.out.dada2_out
            .map { meta, maps, asv_seqs, filt_reads -> [meta, maps] }
            .join(count_taxonomy, failOnDuplicate: true, failOnMismatch: true)
        MAKE_ASV_COUNT_TABLES(map_count_input)
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions)

        if (params.COI_bold_ref_db) {
            LCA_KRONA_REPORTS_BOLD(
                MAKE_ASV_COUNT_TABLES.out.asv_read_counts_out,
                VSEARCH_ASV_LCA_BOLD.out.lca_all,
                VSEARCH_ASV_LCA_BOLD.out.lca_top
            )
            ch_versions = ch_versions.mix(LCA_KRONA_REPORTS_BOLD.out.versions)
        }
        if (params.COI_midori2_ref_db) {
            LCA_KRONA_REPORTS_MIDORI2(
                MAKE_ASV_COUNT_TABLES.out.asv_read_counts_out,
                VSEARCH_ASV_LCA_MIDORI2.out.lca_all,
                VSEARCH_ASV_LCA_MIDORI2.out.lca_top
            )
            ch_versions = ch_versions.mix(LCA_KRONA_REPORTS_MIDORI2.out.versions)
        }
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CLEAN.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(READS_QC.out.fastp_summary_json.map { it[1] })
    ch_versions = ch_versions.mix(FASTQC_CLEAN.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'envident_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    /*****************************/
    /* End of execution reports */
    /****************************/
 
    // Extract runs that failed SeqFu check //
    READS_QC.out.seqfu_check
        .splitCsv(sep: "\t", elem: 1)
        .filter { meta, seqfu_res ->
            seqfu_res[0] != "OK"
        }
        .map { meta, _seqfu_res -> "${meta.id},seqfu_fail" }
        .set { seqfu_fails }

    // Extract runs that failed Suffix Header check //
    READS_QC.out.suffix_header_check
        .filter { meta, sfxhd_res ->
            sfxhd_res.countLines() != 0
        }
        .map { meta, _sfxhd_res -> "${meta.id},sfxhd_fail"  }
        .set { sfxhd_fails }

    // Extract runs that failed Library Strategy check //
    READS_QC.out.amplicon_check
        .filter { meta, strategy ->
            strategy != "AMPLICON"
        }
        .map { meta, _strategy -> "${meta.id},libstrat_fail" }
        .set { libstrat_fails }

    // Extract runs that had zero reads after fastp //
    extended_reads_qc.qc_fail.map { meta, _no_reads -> "${meta.id},min_reads"  }
        .set { min_reads_fails }

    // Extract runs that failed the reads_percentage_threshold parameter
    PROFILE_HMMSEARCH_PFAM.out.profile
        .map { meta, tsv_file -> [meta.id, meta] }
        .join(ch_passed_samples.map { meta -> [meta.id, meta] }, remainder: true)
        .filter { id, meta, passed_meta ->
            passed_meta == null  // These are samples that didn't pass
        }
        .map { id, meta, _low_percent -> "${meta.id},reads_percentage_fail" }
        .set { reads_percentage_fails }

    reads_before_hmm.qc_empty
        .map { meta, _reads -> "${meta.id},empty_after_qc" }
        .set { empty_after_qc_fails }

    // Classify DADA2 results, giving the failure flag precedence over available reports.
    DADA2_SWF.out.dada2_report.map { meta, _dada2_report -> [["id": meta.id, "single_end": meta.single_end], "has_dada2_report"] }
    .concat(
        ch_passed_samples.map { meta -> [["id": meta.id, "single_end": meta.single_end], "qc_pass"] },
        dada2_stats_fail
    )
    .groupTuple()
    .branch { meta, results ->
        failed: results.contains("true")
        passed: results.contains("qc_pass") &&
            results.contains("has_dada2_report") &&
            results.contains("false")
        no_asvs: true
    }
    .set { dada2_qc_results }

    dada2_qc_results.failed
        .map { meta, _results -> "${meta.id},dada2_stats_fail" }
        .set { dada2_failed_runs }

    dada2_qc_results.no_asvs
        .map { meta, _results -> "${meta.id},no_asvs" }
        .set { no_asvs_failed_runs }

    // Save all failed runs to file //
    all_failed_runs = seqfu_fails.concat( sfxhd_fails, libstrat_fails, min_reads_fails, reads_percentage_fails, empty_after_qc_fails, dada2_failed_runs, no_asvs_failed_runs)
    all_failed_runs.collectFile(name: "qc_failed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

    dada2_qc_results.passed
        .map { meta, _results -> "${meta.id},all_results" }
        .set { final_passed_runs }

    // Save all passed runs to file //
    final_passed_runs.collectFile(name: "qc_passed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)
    .set { passed_runs_path }


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
