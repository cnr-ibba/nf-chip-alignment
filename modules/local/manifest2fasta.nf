
process MANIFEST2FASTA {
    tag "$manifest"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biopython=1.78" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    path(manifest)

    output:
    path("*.fasta"), emit: fasta

    script:
    def args = task.ext.args ?: ''

    """
    chip2fasta.py \\
        $args \\
        --input $manifest \\
        --output ${manifest.baseName}.fasta
    """
}
