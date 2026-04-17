include { FASTAEMBEDLENGTH } from '../../../modules/local/fastaembedlength/main'
include { SEQKIT_TRANSLATE } from '../../../modules/nf-core/seqkit/translate/main'
include { HMMER_HMMSEARCH } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { PARSEHMMSEARCHCOVERAGE } from '../../../modules/local/parsehmmsearchcoverage/main'
include { COMBINEHMMSEARCHTBL } from '../../../modules/local/combinehmmsearchtbl/main'

workflow PROFILE_HMMSEARCH_PFAM {

    take:
    reads_fasta
    pfam_db
    reads_json

    main:
    ch_versions = Channel.empty()

    read_counts = reads_json.map { meta, json_file ->
        def json = new groovy.json.JsonSlurper().parseText(json_file.text)
        tuple(meta, json["summary"]["after_filtering"]["total_reads"])
    }
    
    fasta_with_counts = reads_fasta
        .join(read_counts, by: 0)
        .map { meta, fasta, read_count ->
            tuple(meta + ['read_count': read_count], fasta)
        }

    FASTAEMBEDLENGTH(fasta_with_counts, file("${projectDir}/bin/fastx_embed_length.py"))

    SEQKIT_TRANSLATE(FASTAEMBEDLENGTH.out.fasta)

    /*
    * Adaptive FASTA chunking for HMMER parallelization
    * 
    * Uses ideal chunk size (better performance) for small files,
    * but switches to larger chunks for big files to respect the
    * maximum chunk limit (prevents job scheduler overload).
    * 
    * Sequence-based splitting ensures predictable processing times
    * and prevents truncated sequences.
    */

    ch_chunked_pfam_in = SEQKIT_TRANSLATE.out.fastx
        .flatMap{ meta, fasta ->
            def totalSeqs = fasta.countFasta()
            def maxChunks = params.hmmsearch_max_chunks ?: 200
            def idealSeqsPerChunk = params.hmmsearch_seqs_per_chunk ?: 100000

            // Adaptive strategy: ideal chunks for small files, larger chunks for big files
            def seqsPerChunk = (totalSeqs <= maxChunks * idealSeqsPerChunk) ? 
                idealSeqsPerChunk :                           
                Math.ceil(totalSeqs / maxChunks) as Integer   
                
            def chunks = fasta.splitFasta(by: seqsPerChunk, file: true)
            
            chunks.collect{ chunk -> tuple(groupKey(meta, chunks.size()), chunk) }
        }
        .combine(pfam_db)
        .map{ meta, reads, db -> [meta, db, reads, false, true, true] }

    HMMER_HMMSEARCH(ch_chunked_pfam_in)
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    COMBINEHMMSEARCHTBL(
        HMMER_HMMSEARCH.out.domain_summary.groupTuple()
    )

    PARSEHMMSEARCHCOVERAGE(COMBINEHMMSEARCHTBL.out.concatenated_result, file("${projectDir}/bin/hmmer_domtbl_parse_coverage.py"))

    def combined_hmm = COMBINEHMMSEARCHTBL.out.concatenated_result.map { meta, reads ->
        [[id: meta.id, single_end: meta.single_end], reads]
    }

    def coverage_tsv = PARSEHMMSEARCHCOVERAGE.out.tsv.map { meta, reads ->
        [[id: meta.id, single_end: meta.single_end], reads]
    }

    ch_versions = ch_versions.mix(PARSEHMMSEARCHCOVERAGE.out.versions)

    emit:
    profile  = coverage_tsv
    domtbl   = combined_hmm
    versions = ch_versions                                   // channel: [ versions.yml ]
}
