
process CUTADAPT {
    tag "$meta.id"
    label 'very_light'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py311haab0aaa_0':
        'biocontainers/cutadapt:5.2--py311haab0aaa_0' }"

    input:
    tuple val(meta), path(reads), path(primers)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    tuple val(meta), path('*.json')         , emit: json
    tuple val(meta), path('*_summ.tsv')     , optional: true, emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def first_trimmed = meta.single_end
        ? "-o ${prefix}.first_pass.fastq.gz --info-file ${prefix}.first_pass.tsv"
        : "-o ${prefix}_1.first_pass.fastq.gz -p ${prefix}_2.first_pass.fastq.gz --info-file ${prefix}_1.first_pass.tsv --info-file-paired ${prefix}_2.first_pass.tsv"
    def final_trimmed = meta.single_end
        ? "-o ${prefix}.trim.fastq.gz --info-file ${prefix}.final.tsv"
        : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz --info-file ${prefix}_1.final.tsv --info-file-paired ${prefix}_2.final.tsv"
    def first_primer_args = primers[0].size() > 0 ? "-g file:${primers[0]}" : ""
    if (!meta.single_end && primers[1].size() > 0) {
        first_primer_args += " -G file:${primers[1]}"
    }
    def second_primer_args = primers[3].size() > 0 ? "-a file:${primers[3]}" : ""
    if (!meta.single_end && primers[2].size() > 0) {
        second_primer_args += " -A file:${primers[2]}"
    }

    // Produces summary tsvs of the positions primers were trimmed at
    def stages = ['first_pass', 'final']

    def stats_cmd = stages.collect { stage ->

        def input_files = meta.single_end
            ? ["${prefix}.${stage}.tsv"]
            : ["${prefix}_1.${stage}.tsv", "${prefix}_2.${stage}.tsv"]

        def output_files = meta.single_end
            ? ["${prefix}_${stage}_se_summ.tsv"]
            : ["${prefix}_${stage}_1_summ.tsv", "${prefix}_${stage}_2_summ.tsv"]

        input_files
            .withIndex()
            .collect { input_file, i ->
                """
                awk -F'\\t' '{print length(\$5),length(\$6),length(\$7)}' ${input_file} |
                    sort |
                    uniq -c |
                    awk 'OFS="\\t" {print \$2,\$3,\$4,\$1}' |
                    sort -nr -k4,4 > ${output_files[i]}
                """
            }
            .join('\n')
    }
    .join('\n')

    if(first_primer_args == ""){
        if (!meta.single_end){
            """
            touch ${prefix}.cutadapt.log
            touch ${prefix}_1.trim.fastq.gz
            touch ${prefix}_2.trim.fastq.gz
            touch ${prefix}.cutadapt.json

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(cutadapt --version)
            END_VERSIONS
            """
        }
        else{
            """
            touch ${prefix}.cutadapt.log
            touch ${prefix}.trim.fastq.gz
            touch ${prefix}.cutadapt.json

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cutadapt: \$(cutadapt --version)
            END_VERSIONS
            """

        }
    }
    else{
        """
        cutadapt \\
            --cores $task.cpus \\
            $args \\
            ${first_trimmed} \\
            ${first_primer_args} \\
            ${reads} \\
            --json ${prefix}.first_pass.cutadapt.json \\
            > ${prefix}.first_pass.cutadapt.log

        cutadapt \\
            --cores $task.cpus \\
            $args \\
            ${final_trimmed} \\
            ${second_primer_args} \\
            ${meta.single_end ? "${prefix}.first_pass.fastq.gz" : "${prefix}_1.first_pass.fastq.gz ${prefix}_2.first_pass.fastq.gz"} \\
            --json ${prefix}.cutadapt.json \\
            > ${prefix}.cutadapt.log
        
        ${stats_cmd}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end
        ? "${prefix}.trim.fastq.gz"
        : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    def summaries = meta.single_end
        ? "${prefix}_first_pass_se_summ.tsv ${prefix}_final_se_summ.tsv"
        : "${prefix}_first_pass_1_summ.tsv ${prefix}_first_pass_2_summ.tsv ${prefix}_final_1_summ.tsv ${prefix}_final_2_summ.tsv"
    """
    touch ${prefix}.cutadapt.log
    touch ${prefix}.cutadapt.json
    touch ${trimmed}
    touch ${summaries}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
