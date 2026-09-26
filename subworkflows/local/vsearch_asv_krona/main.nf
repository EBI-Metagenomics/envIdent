include { VSEARCH_ASV_LCA } from '../vsearch_asv_lca/main'
include { LCA_KRONA_REPORTS } from '../lca_krona_reports/main'

workflow VSEARCH_ASV_KRONA {
    take:
        dada2_output
        ref_db
        asv_counts // Shared count table generated once after DADA2

    main:
        VSEARCH_ASV_LCA(dada2_output, ref_db)
        LCA_KRONA_REPORTS(
            asv_counts,
            VSEARCH_ASV_LCA.out.lca_all,
            VSEARCH_ASV_LCA.out.lca_top
        )

    emit:
        vsearch_out = VSEARCH_ASV_LCA.out.vsearch_out
        lca_input = VSEARCH_ASV_LCA.out.lca_input
        clean_hits = VSEARCH_ASV_LCA.out.clean_hits
        lca_all = LCA_KRONA_REPORTS.out.lca_all
        lca_top = LCA_KRONA_REPORTS.out.lca_top
        krona_all_counts = LCA_KRONA_REPORTS.out.krona_all_counts
        krona_top_counts = LCA_KRONA_REPORTS.out.krona_top_counts
        krona_all = LCA_KRONA_REPORTS.out.krona_all
        krona_top = LCA_KRONA_REPORTS.out.krona_top
        versions = VSEARCH_ASV_LCA.out.versions.mix(LCA_KRONA_REPORTS.out.versions)
}
