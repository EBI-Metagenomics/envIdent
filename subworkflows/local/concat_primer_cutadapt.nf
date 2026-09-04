
include { REV_COMP_SE_PRIMERS     } from '../../modules/local/rev_comp_se_primers/main.nf'
include { CUTADAPT                } from '../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        reads
    main:

        ch_versions = Channel.empty()
        
        REV_COMP_SE_PRIMERS(
            concat_input
        )
        ch_versions = ch_versions.mix(REV_COMP_SE_PRIMERS.out.versions.first())

        // Join separate forward and reverse primers to the fastp-cleaned reads and run cutadapt.

        cutadapt_input = REV_COMP_SE_PRIMERS.out.rev_comp_se_primers_out
            .map { meta, fwd_primer, rev_primer ->
                [meta.subMap('id', 'single_end'), meta, [fwd_primer, rev_primer]]
            }
            .join(
                reads.map { meta, final_reads ->
                    [meta.subMap('id', 'single_end'), final_reads]
                },
                by: [0]
            )
            .map { _key, meta, primers, final_reads ->
                [meta, final_reads, primers]
            }

        CUTADAPT(
            cutadapt_input
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    emit:
        cutadapt_out = CUTADAPT.out.reads
        cutadapt_json = CUTADAPT.out.json
        versions = ch_versions

}