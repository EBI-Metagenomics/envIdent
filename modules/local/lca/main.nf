process LCA {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://docker.io/library/node:24'
        : 'docker.io/library/node:24' }"

    input:
    tuple val(meta), path(vsearch_lca)

    output:
    tuple val(meta), path("*_lca_all.tsv"), emit: lca_all
    tuple val(meta), path("*_lca_top.tsv"), emit: lca_top    
    path "versions.yml", emit: versions

    script:
    """
    all_lowest_common_ancestor.js ${vsearch_lca} > ${meta.id}_lca_all.tsv
    top_lowest_common_ancestor.js ${vsearch_lca} > ${meta.id}_lca_top.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        node: \$(node --version |& sed 's/v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_lca_all.tsv
    touch ${meta.id}_lca_top.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        node: \$(node --version |& sed 's/v//')
    END_VERSIONS
    """
}
