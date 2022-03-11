#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'

// get manifest from parameters
manifest_ch = Channel.fromPath( params.manifest )


process MANIFEST2FASTA {
    tag 'fastaconversion'
    label 'process_low'

    container = 'images/chip-conversion.sif'

    input:
    path(manifest)

    output:
    path("*.fasta"), emit: fasta

    """
    chip2fasta.py --input $manifest --output ${manifest.baseName}.fasta
    """
}

workflow {
    MANIFEST2FASTA(manifest_ch).view()
}
