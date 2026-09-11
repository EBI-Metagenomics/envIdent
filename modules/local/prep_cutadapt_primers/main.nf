process PREP_CUTADAPT_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.84'
        : 'quay.io/biocontainers/biopython:1.84'}"

    input:
    tuple val(meta), path(forward_primers), path(reverse_primers)

    output:
    tuple val(meta),
        path("${meta.id}_F.fasta"),
        path("${meta.id}_R.fasta"),
        path("${meta.id}_F_RC.fasta"),
        path("${meta.id}_R_RC.fasta"),
        emit: cutadapt_prepared_primers
    path "versions.yml", emit: versions

    script:
    def forward_rc = "${meta.id}_F_rc_notprep.fasta"
    def reverse_rc = "${meta.id}_R_rc_notprep.fasta"
    def forward_prepared = "${meta.id}_F.fasta"
    def reverse_prepared = "${meta.id}_R.fasta"
    def forward_rc_prepared = "${meta.id}_F_RC.fasta"
    def reverse_rc_prepared = "${meta.id}_R_RC.fasta"

    """
    
    revcomp_primers.py ${forward_primers} ${forward_rc}
    revcomp_primers.py ${reverse_primers} ${reverse_rc}

    prep_cutadapt_primers.py --position start --replace-i --n ${params.bases_allowed_before_primer} ${forward_primers} ${forward_prepared}
    prep_cutadapt_primers.py --position start --replace-i --n ${params.bases_allowed_before_primer} ${reverse_primers} ${reverse_prepared}
    prep_cutadapt_primers.py --position end --replace-i --n ${params.bases_allowed_before_primer} ${forward_rc} ${forward_rc_prepared}
    prep_cutadapt_primers.py --position end --replace-i --n ${params.bases_allowed_before_primer} ${reverse_rc} ${reverse_rc_prepared}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

        stub:
    """
    touch ${meta.id}_F.fasta
    touch ${meta.id}_R.fasta
    touch ${meta.id}_F_RC.fasta
    touch ${meta.id}_R_RC.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
