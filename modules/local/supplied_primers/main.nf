process SUPPLIED_PRIMERS {
    tag "$meta.id"
    label 'very_light'

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}_provided_primers.fasta"), emit: supplied_primer_out

    script:
    """
    cat > ${meta.id}_provided_primers.fasta <<-END_PRIMERS
    >${meta.id}_provided|F
    ${meta.forward_primer}
    >${meta.id}_provided|R
    ${meta.reverse_primer}
    END_PRIMERS
    """
}