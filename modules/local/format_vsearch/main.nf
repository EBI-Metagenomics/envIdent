process FORMAT_VSEARCH {
    tag "$meta.id"
    label 'very_light'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.84'
        : 'quay.io/biocontainers/biopython:1.84'}"

    input:
    tuple val(meta), path(vsearch_out)
    val mode 

    output:
    tuple val(meta), path("${meta.id}_${mode}_formatted.tsv"), emit: formatted    
    path "versions.yml", emit: versions

    script:
    def outfile = "${meta.id}_${mode}_formatted.tsv"
    def format

    if (mode == 'lca') {
        format = "formatting_vsearch.py -i ${vsearch_out} -o ${outfile} -p -q --lca"
    } else if (mode == 'clean') {
        format = "formatting_vsearch.py -i ${vsearch_out} -o ${outfile} -a -p"
    } else {
        error "Unknown FORMAT_VSEARCH mode '${mode}': expected 'lca' or 'clean'"
    }

    """
    ${format}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${mode}_formatted.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
