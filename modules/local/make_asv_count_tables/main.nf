
process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(maps), path(asvtaxtable), path(reads)
    val db_label

    output:
    tuple val(meta), path("*_asv_read_counts.tsv"), optional: true, emit: asv_read_counts_out
    tuple val(meta), path("*_asv_krona_counts.txt"), optional: true, emit: asv_count_tables_out
    path "versions.yml"                           , emit: versions

    script:
    """
    if [[ ${meta.single_end} = true ]]; then
        zcat ${reads} | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table.py -t ${asvtaxtable} -f ${maps} -hd headers.txt -s ${meta.id}_${meta.var_region}_${db_label}
    else
        zcat ${reads[0]} | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table.py -t ${asvtaxtable} -f ${maps[0]} -r ${maps[1]} -hd headers.txt -s ${meta.id}_${meta.var_region}_${db_label}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.var_region}_${db_label}_asv_read_counts.tsv
    touch ${meta.id}_${meta.var_region}_${db_label}_asv_krona_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
