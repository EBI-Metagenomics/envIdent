/**
 * Join QC and trimmed reads by id/single_end, returning [meta, reads] tuples.
 * Both inputs carry [Map meta, Path reads] for single-end data or
 * [Map meta, List<Path> reads] for paired-end data. Unmatched keys are omitted.
 *
 * Preserve the existing fallback: use QC reads when the inspected trimmed file
 * is zero bytes. This checks file size, not FASTQ record count; an empty gzip
 * stream still has a nonzero size. Only id/single_end metadata is retained.
 */
def reads_merged_input_prep(reads_qc, cutadapt_channel) {

    def selected_reads = reads_qc
        .map { meta, reads -> [ meta.subMap('id', 'single_end'), reads ] }
        .join(
            cutadapt_channel.map { meta, reads -> [ meta.subMap('id', 'single_end'), reads ] },
            by: 0
        )

        .map { meta, fastp_reads, cutadapt_reads ->

            def trimmed_file_size = meta.single_end
                ? cutadapt_reads.size()
                : cutadapt_reads[0].size()

            def final_reads = trimmed_file_size > 0 ? cutadapt_reads : fastp_reads

            [ meta, final_reads ]
    }

    return selected_reads
}
