process PARSESAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0'
        : 'biocontainers/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0'}"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*_parsed.tsv"), emit: parsed
    tuple val(meta), path("*_parsed.sam"), emit: sam
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parseSAM.py \\
        ${args} \\
        -i ${sam} \\
        -o ${prefix}_parsed.tsv \\
        -O ${prefix}_parsed.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def _args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_parsed.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
