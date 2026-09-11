process EXTRACT_CUTADAPT_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0'
        : 'biocontainers/seqkit:2.13.0--he881be0_0'}"

    input:
    tuple val(meta), path(forward_detected), path(reverse_detected)
    path(cutadapt_primer_library, stageAs: 'cutadapt_primer_library')

    output:
    tuple val(meta),
        path("${meta.id}_F.fasta"),
        path("${meta.id}_R.fasta"),
        path("${meta.id}_F_RC.fasta"),
        path("${meta.id}_R_RC.fasta"),
        emit: cutadapt_prepared_primers
    path "versions.yml", emit: versions

    script:
    def forward_prepared = "${meta.id}_F.fasta"
    def reverse_prepared = "${meta.id}_R.fasta"
    def forward_rc_prepared = "${meta.id}_F_RC.fasta"
    def reverse_rc_prepared = "${meta.id}_R_RC.fasta"

    """
    seqkit grep -n -f <(seqkit seq -n ${forward_detected}) ${cutadapt_primer_library}/*F.fasta > ${forward_prepared}
    seqkit grep -n -f <(seqkit seq -n ${reverse_detected}) ${cutadapt_primer_library}/*R.fasta > ${reverse_prepared}
    seqkit grep -n -f <(seqkit seq -n ${forward_detected}) ${cutadapt_primer_library}/*F_RC.fasta > ${forward_rc_prepared}
    seqkit grep -n -f <(seqkit seq -n ${reverse_detected}) ${cutadapt_primer_library}/*R_RC.fasta > ${reverse_rc_prepared}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primer-preparation: python
    END_VERSIONS
    """

    stub:
    def forward_prepared = "${meta.id}_F.fasta"
    def reverse_prepared = "${meta.id}_R.fasta"
    def forward_rc_prepared = "${meta.id}_F_RC.fasta"
    def reverse_rc_prepared = "${meta.id}_R_RC.fasta"

    """
    touch ${forward_prepared}
    touch ${reverse_prepared}
    touch ${forward_rc_prepared}
    touch ${reverse_rc_prepared}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primer-preparation: python
    END_VERSIONS
    """
}