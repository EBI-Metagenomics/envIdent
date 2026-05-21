import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

def reads_merged_input_prep( reads_qc, cutadapt_channel ) {

    def dada2_input = reads_qc
        .map { meta, reads -> [ meta.subMap('id', 'single_end'), reads ] }
        .join(
            cutadapt_channel.map { meta, reads -> [ meta.subMap('id', 'single_end'), reads ] },
            by: 0
        )

        .map { meta, fastp_reads, cutadapt_reads ->

            def cutadapt_read_size = meta.single_end
                ? cutadapt_reads.size()
                : cutadapt_reads[0].size()
        
            def final_reads = cutadapt_read_size > 0 ? cutadapt_reads : fastp_reads
        
            [ meta, final_reads ]
    }

    return dada2_input
}
