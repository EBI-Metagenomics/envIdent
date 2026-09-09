include { CUTADAPT                } from '../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    take:
        concat_input
        reads
    main:
        ch_versions = Channel.empty()

        // Join four prepared primer files to the fastp-cleaned reads.
        cutadapt_input = concat_input
            .map { meta, fwd_primer, rev_primer, fwd_primer_rc, rev_primer_rc ->
                [meta.subMap('id', 'single_end'), meta, [fwd_primer, rev_primer, fwd_primer_rc, rev_primer_rc]]
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
