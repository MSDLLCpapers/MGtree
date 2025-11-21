process MGTREE_COUNT {
    tag "${meta.id}"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0'
        : 'biocontainers/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0'}"

    input:
    tuple val(meta), path(tsv)
    path newick
    val no_leaf_genotypes

    output:
    tuple val(meta), path("*_genotype.tsv"), emit: genotype, optional: true
    tuple val(meta), path("*_tally.tsv"), emit: tally, optional: true
    tuple val(meta), stdout, emit: console
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def g_flag = no_leaf_genotypes ? '-g' : ''
    """
    MGtree.py \\
        ${args} \\
        ${g_flag} \\
        -i ${tsv} \\
        -o ${prefix}_genotype.tsv \\
        -t ${prefix}_tally.tsv \\
        -n ${newick}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """

    stub:
    def _args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_genotype.tsv \\
          ${prefix}_tally.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
        matplotlib: \$(python -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """
}
