process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(maps), path(reads)

    output:
    tuple val(meta), path('*_asv_read_counts.tsv'), emit: asv_read_counts
    path 'versions.yml', emit: versions

    script:
    def forward_map = meta.single_end ? maps : maps[0]
    def forward_reads = meta.single_end ? reads : reads[0]
    """
    # DADA2 input has already been restricted to marker-matching reads.
    # Use every filtered-read ID as the toolkit's allowed-read list.
    # The toolkit splits headers on a space; include one for bare read IDs too.
    zcat ${forward_reads} | awk 'NR % 4 == 1 {print \$1 " "}' > headers.txt
    awk '{sub(/^@/, "", \$1); print \$1}' headers.txt > reads.COI.txt
    printf 'ASV\\tSuperkingdom\\tKingdom\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\n' > empty_tax.tsv
    make_asv_count_table -t empty_tax.tsv -f ${forward_map} -a reads.COI.txt -hd headers.txt -s ${meta.id}
    mv ${meta.id}_COI.txt_asv_read_counts.tsv ${meta.id}_asv_read_counts.tsv
    if [[ ! -s ${meta.id}_asv_read_counts.tsv ]]; then
        printf 'asv\\tcount\\n' >> ${meta.id}_asv_read_counts.tsv
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

    stub:
    """
    printf 'asv\\tcount\\n' > ${meta.id}_asv_read_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
