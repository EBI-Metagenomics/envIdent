include { VSEARCH_USEARCHGLOBAL as VSEARCH       } from '../../../modules/nf-core/vsearch/usearchglobal/main'
include { FORMAT_VSEARCH as FORMAT_VSEARCH_LCA   } from '../../../modules/local/format_vsearch/main'
include { FORMAT_VSEARCH as FORMAT_VSEARCH_CLEAN } from '../../../modules/local/format_vsearch/main'
include { LCA                                    } from '../../../modules/local/lca/main'

workflow VSEARCH_ASV_LCA {
    
    take:
        vsearch_input // [meta, asv_seqs]
        ref_db        // Reference FASTA with taxonomy in sequence headers

    main:

        ch_versions = channel.empty()

        VSEARCH(
            vsearch_input,
            ref_db,
            params.idcutoff,
            'userout',
            "query+target+id+qcov+alnlen"
        )
        ch_versions = ch_versions.mix(VSEARCH.out.versions.first())

        FORMAT_VSEARCH_LCA(
            VSEARCH.out.tsv,
            'lca'
        )
        ch_versions = ch_versions.mix(FORMAT_VSEARCH_LCA.out.versions.first())

        FORMAT_VSEARCH_CLEAN(
            VSEARCH.out.tsv,
            'clean'
        )
        ch_versions = ch_versions.mix(FORMAT_VSEARCH_CLEAN.out.versions.first())

        LCA(FORMAT_VSEARCH_LCA.out.formatted)
        ch_versions = ch_versions.mix(LCA.out.versions.first())

    emit:
        vsearch_out = VSEARCH.out.tsv
        lca_input = FORMAT_VSEARCH_LCA.out.formatted
        clean_hits = FORMAT_VSEARCH_CLEAN.out.formatted
        lca_all = LCA.out.lca_all
        lca_top = LCA.out.lca_top
        versions = ch_versions
}
