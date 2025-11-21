process NGMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ngmerge:0.3--ha92aebf_1'
        : 'biocontainers/ngmerge:0.3--ha92aebf_1'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fq.gz"), emit: merged_reads, optional: true
    tuple val(meta), path("*_1.fastq.gz"), emit: unstitched_read1
    tuple val(meta), path("*_2.fastq.gz"), emit: unstitched_read2
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        error("NGmerge is not compatible with single-end samples")
    }

    def (r1, r2) = reads.collate(2).transpose().collect { r -> r.join(',') }

    """
    NGmerge \\
        -1 ${r1} \\
        -2 ${r2} \\
        -o ${prefix}.merged.fq.gz \\
        -f ${prefix}_unstitched \\
        -z \\
        -n ${task.cpus} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$(echo \$(NGmerge --version 2>&1) | sed 's/^.*NGmerge, version //; s/ Copyright.*// ; s/: //g' ))
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.merged.fq.gz
    touch ${prefix}.merged.fq.gz_1.fastq.gz
    touch ${prefix}.merged.fq.gz_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGmerge: \$(echo \$(NGmerge --version 2>&1) | sed 's/^.*NGmerge, version //; s/ Copyright.*// ; s/: //g' ))
    END_VERSIONS
    """
}
