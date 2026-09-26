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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // primers[] convention: [forward, reverse, forward_RC, reverse_RC]
    def (fwd_primer, rev_primer, fwd_primer_rc, rev_primer_rc) = primers

    def first_primer_args  = fwd_primer.size() > 0 ? "-g file:${fwd_primer}" : ''
    if (!meta.single_end && rev_primer.size() > 0) {
        first_primer_args += " -G file:${rev_primer}"
    }
    def second_primer_args = rev_primer_rc.size() > 0 ? "-a file:${rev_primer_rc}" : ''
    if (!meta.single_end && fwd_primer_rc.size() > 0) {
        second_primer_args += " -A file:${fwd_primer_rc}"
    }

    // Two sequential passes: primers are trimmed and independently logged/
    // stat'd for the 5' pass and the 3' pass. cutadapt is happy to receive
    // an empty primer arg for a pass that has nothing to trim, so we only
    // skip cutadapt entirely when NEITHER pass has any primers at all.
    // Each pass carries everything needed to build its own cutadapt call,
    // stats step, and "no primers" placeholder outputs, so nothing below
    // has to be retyped per pass.
    def passes = [
        [
            stage      : 'first_pass',
            primer_args: first_primer_args,
            input      : "${reads}",
            out_files  : meta.single_end
                ? ["${prefix}.first_pass.fastq.gz"]
                : ["${prefix}_1.first_pass.fastq.gz", "${prefix}_2.first_pass.fastq.gz"],
            tsv_files  : meta.single_end
                ? ["${prefix}.first_pass.tsv"]
                : ["${prefix}_1.first_pass.tsv", "${prefix}_2.first_pass.tsv"],
            summ_files : meta.single_end
                ? ["${prefix}_first_pass_se_summ.tsv"]
                : ["${prefix}_first_pass_1_summ.tsv", "${prefix}_first_pass_2_summ.tsv"],
        ],
        [
            stage      : 'second_pass',
            primer_args: second_primer_args,
            input      : meta.single_end
                ? "${prefix}.first_pass.fastq.gz"
                : "${prefix}_1.first_pass.fastq.gz ${prefix}_2.first_pass.fastq.gz",
            out_files  : meta.single_end
                ? ["${prefix}.trim.fastq.gz"]
                : ["${prefix}_1.trim.fastq.gz", "${prefix}_2.trim.fastq.gz"],
            tsv_files  : meta.single_end
                ? ["${prefix}.second_pass.tsv"]
                : ["${prefix}_1.second_pass.tsv", "${prefix}_2.second_pass.tsv"],
            summ_files : meta.single_end
                ? ["${prefix}_second_pass_se_summ.tsv"]
                : ["${prefix}_second_pass_1_summ.tsv", "${prefix}_second_pass_2_summ.tsv"],
        ],
    ]

    def no_primers_at_all = (first_primer_args == '' && second_primer_args == '')

    if (no_primers_at_all) {
        """
        ${passes.collect { p -> """
        touch ${prefix}.${p.stage}.cutadapt.log
        touch ${p.out_files.join(' ')}
        echo '{"message": "No primers were inputted so trimming not performed"}' > ${prefix}.${p.stage}.cutadapt.json
        ${p.summ_files.collect { f -> "printf 'Pre-primer length\\tPrimer length\\tPost-primer length\\tRead count\\n' > ${f}" }.join('\n')}
        """ }.join('\n')}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        ${passes.collect { p ->
            def io = meta.single_end
                ? "-o ${p.out_files[0]} --info-file ${p.tsv_files[0]}"
                : "-o ${p.out_files[0]} -p ${p.out_files[1]} --info-file ${p.tsv_files[0]} --info-file-paired ${p.tsv_files[1]}"
            """
        cutadapt \\
            --cores $task.cpus \\
            $args \\
            ${io} \\
            ${p.primer_args} \\
            ${p.input} \\
            --json ${prefix}.${p.stage}.cutadapt.json \\
            > ${prefix}.${p.stage}.cutadapt.log
        """
        }.join('\n')}

        ${passes.collect { p ->
            p.tsv_files.withIndex().collect { tsv, i -> "primer_trim_stats.sh ${tsv} ${p.summ_files[i]}" }.join('\n')
        }.join('\n')}

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
        ? "${prefix}_first_pass_se_summ.tsv ${prefix}_second_pass_se_summ.tsv"
        : "${prefix}_first_pass_1_summ.tsv ${prefix}_first_pass_2_summ.tsv ${prefix}_second_pass_1_summ.tsv ${prefix}_second_pass_2_summ.tsv"
    """
    touch ${prefix}.first_pass.cutadapt.log
    touch ${prefix}.second_pass.cutadapt.log
    touch ${prefix}.first_pass.cutadapt.json
    touch ${prefix}.second_pass.cutadapt.json
    touch ${trimmed}
    for summary in ${summaries}; do
        printf 'Pre-primer length\\tPrimer length\\tPost-primer length\\tRead count\\n' > \$summary
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}