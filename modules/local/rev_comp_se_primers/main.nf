
process REV_COMP_SE_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(forward_primers), path(reverse_primers)

    output:
    tuple val(meta), path("${meta.id}_forward.fasta"), path("${meta.id}_reverse.fasta"), emit: rev_comp_se_primers_out
    path "versions.yml"                                , emit: versions

    script:
    """
    if [[ ${meta.single_end} = true ]]; then
        rev_comp_se_primers -i $reverse_primers -s ${meta.id} -o ./
        mv ./${meta.id}_rev_comp_se_primers.fasta ./${meta.id}_reverse.fasta
    else
        cp $reverse_primers ./${meta.id}_reverse.fasta
    fi

    cp $forward_primers ./${meta.id}_forward.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
