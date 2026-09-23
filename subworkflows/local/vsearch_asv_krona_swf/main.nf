
include { VSEARCH_USEARCHGLOBAL as VSEARCH } from '../../../modules/nf-core/vsearch/usearchglobal/main'
include { FORMAT_VSEARCH as FORMAT_VSEARCH_LCA } from '../../../modules/local/format_vsearch/main'
include { FORMAT_VSEARCH as FORMAT_VSEARCH_CLEAN } from '../../../modules/local/format_vsearch/main'

include { LCA } from '../../../modules/local/lca/main'
include { LCA_KRONA_COUNTS } from '../../../modules/local/lca_krona_counts/main'
include { KRONA_KTIMPORTTEXT as KRONA_ALL } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { KRONA_KTIMPORTTEXT as KRONA_TOP } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'

workflow VSEARCH_ASV_LCA {
    
    take:
        dada2_output // [meta, maps, asv_seqs, filtered_reads]
        ref_db       // Reference FASTA with taxonomy in sequence headers

    main:

        ch_versions = channel.empty()

        vsearch_input = dada2_output
                       .map { meta, maps, asv_seqs, filt_reads ->
                            [ meta, asv_seqs ]
                        }

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

workflow LCA_KRONA_REPORTS {
    take:
        asv_counts
        lca_all
        lca_top

    main:
        ch_versions = channel.empty()

        // Join on the complete metadata map to keep samples and marker regions separate.
        counts_input = asv_counts
            .join(lca_all, failOnDuplicate: true, failOnMismatch: true)
            .join(lca_top, failOnDuplicate: true, failOnMismatch: true)
        LCA_KRONA_COUNTS(counts_input)
        ch_versions = ch_versions.mix(LCA_KRONA_COUNTS.out.versions.first())

        KRONA_ALL(LCA_KRONA_COUNTS.out.all_counts)
        KRONA_TOP(LCA_KRONA_COUNTS.out.top_counts)
        ch_versions = ch_versions.mix(KRONA_ALL.out.versions.first(), KRONA_TOP.out.versions.first())

    emit:
        lca_all = LCA_KRONA_COUNTS.out.lca_all
        lca_top = LCA_KRONA_COUNTS.out.lca_top
        krona_all_counts = LCA_KRONA_COUNTS.out.all_counts
        krona_top_counts = LCA_KRONA_COUNTS.out.top_counts
        krona_all = KRONA_ALL.out.html
        krona_top = KRONA_TOP.out.html
        versions = ch_versions
}
