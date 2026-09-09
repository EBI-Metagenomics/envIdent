
process STD_PRIMER_FLAG {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)
    tag "$meta.id"
    label 'med_cpu_light_mem'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mi-pimento:${params.pimento_version}":
        "biocontainers/mi-pimento:${params.pimento_version}"}"
    input:
    tuple val(meta), path(reads)
    path(std_primer_library)

    output:
    tuple val(meta), path("*std_primers.fasta"), emit: std_primer_out
    path "*std_primer_out.txt"
    path "versions.yml"                        , emit: versions

    script:
    def primer_direction = meta.direction.toLowerCase()
    if (!['f', 'r'].contains(primer_direction)) {
        error "Unsupported primer direction '${meta.direction}'; expected 'f' or 'r'"
    }
    def primer_pattern = primer_direction == 'f' ? '*[Ff].fasta' : '*[Rr].fasta'
    def std_primer_library_arg = "${std_primer_library}" ? "-p ${std_primer_library}/${primer_pattern}" : ""
    
    """
    pimento std -i ${reads} ${std_primer_library_arg} --threads $task.cpus -o ${meta.id}_${meta.var_region}_${meta.direction}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """

}
