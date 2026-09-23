process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/pandas:2.2.1":
        "quay.io/biocontainers/pandas:2.2.1" }"

    input:
    tuple val(meta), path(maps), path(asvtaxtable)

    output:
    tuple val(meta), path('*_asv_read_counts.tsv'), emit: asv_read_counts_out
    path 'versions.yml', emit: versions

    script:
    def forward_map = maps instanceof List ? maps.sort { it.name }[0] : maps
    """
    make_asv_count_table.py -t ${asvtaxtable} -f ${forward_map} -s ${meta.id}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python -c 'import platform; print(platform.python_version())')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """

    stub:
    """
    printf 'asv\\tcount\\n' > ${meta.id}_asv_read_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python -c 'import platform; print(platform.python_version())')
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
