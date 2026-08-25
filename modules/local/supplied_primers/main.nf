process SUPPLIED_PRIMERS {
    tag "$meta.id"
    label 'very_light'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}_provided_forward.fasta"), path("${meta.id}_provided_reverse.fasta"), emit: supplied_primer_out

    script:
    """
    cat > ${meta.id}_provided_forward.fasta <<-END_PRIMERS
    >${meta.id}_provided|F
    ${meta.forward_primer}
    END_PRIMERS

    cat > ${meta.id}_provided_reverse.fasta <<-END_PRIMERS
    >${meta.id}_provided|R
    ${meta.reverse_primer}
    END_PRIMERS
    """
}