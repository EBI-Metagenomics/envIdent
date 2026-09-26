process LCA_KRONA_COUNTS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://docker.io/library/node:24'
        : 'docker.io/library/node:24' }"
        
    input:
    tuple val(meta), path(counts), path(lca_all), path(lca_top)

    output:
    tuple val(meta), path('*_lca_all_krona_counts.txt'), emit: all_counts
    tuple val(meta), path('*_lca_top_krona_counts.txt'), emit: top_counts
    tuple val(meta), path('*_lca_all_with_counts.tsv'), emit: lca_all
    tuple val(meta), path('*_lca_top_with_counts.tsv'), emit: lca_top
    path 'versions.yml', emit: versions

    script:
    """
    lca_krona_counts.js ${counts} ${lca_all} ${meta.id}_lca_all_krona_counts.txt ${meta.id}_lca_all_with_counts.tsv
    lca_krona_counts.js ${counts} ${lca_top} ${meta.id}_lca_top_krona_counts.txt ${meta.id}_lca_top_with_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        node: \$(node --version | sed 's/^v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_lca_all_krona_counts.txt ${meta.id}_lca_top_krona_counts.txt ${meta.id}_lca_all_with_counts.tsv ${meta.id}_lca_top_with_counts.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        node: \$(node --version | sed 's/^v//')
    END_VERSIONS
    """
}
