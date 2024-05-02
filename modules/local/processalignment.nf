
process PROCESSALIGNMENT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(txt)

    output:
    path "*.csv", emit: csv
    path "*.aln", emit: aln
    path "*.err", emit: err

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    process_alignment.py \\
        -a ${txt} \\
        -o ${prefix}.csv \\
        $args \\
        --output_aln ${prefix}.aln \\
        --error_csv ${prefix}.err
    """
}
