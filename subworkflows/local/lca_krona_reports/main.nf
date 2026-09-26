include { LCA_KRONA_COUNTS } from '../../../modules/local/lca_krona_counts/main'
include { KRONA_KTIMPORTTEXT as KRONA_ALL } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { KRONA_KTIMPORTTEXT as KRONA_TOP } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'

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
