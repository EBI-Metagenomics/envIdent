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

            // Detect if cutadapt reads exist (SE vs PE check)
            def cutadapt_read_size = meta.single_end
                                        ? cutadapt_reads.size()
                                        : cutadapt_reads[0].size()

            // If cutadapt reads are empty, stage fastp reads into a temp directory
            // to avoid filename collision when the next fastp process runs
            def final_reads
            if (cutadapt_read_size > 0) {
                final_reads = cutadapt_reads
            } else {
                def staging_dir = Files.createTempDirectory("fastp_staging")
                final_reads = fastp_reads instanceof List
                    ? fastp_reads.collect { file ->
                        def path = file instanceof Path ? file : Paths.get(file.toString())
                        def name = path.getFileName().toString().replaceFirst(/\.fastp/, '')
                        def newPath = staging_dir.resolve(name)
                        Files.createSymbolicLink(newPath, path)
                        newPath
                    }
                    : {
                        def path = fastp_reads instanceof Path ? fastp_reads : Paths.get(fastp_reads.toString())
                        def name = path.getFileName().toString().replaceFirst(/\.fastp/, '')
                        def newPath = staging_dir.resolve(name)
                        Files.createSymbolicLink(newPath, path)
                        [ newPath ]
                    }()
            }

            [ meta, final_reads ]
        }

    return dada2_input
}
