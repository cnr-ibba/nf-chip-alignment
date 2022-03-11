#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { UCSC_BLAT } from './modules/local/ucsc/blat'

// get manifest from parameters
manifest_ch = Channel.fromPath( params.manifest )
genome_ch = Channel.fromPath( params.genome )
ch_versions = Channel.empty()

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
    MANIFEST2FASTA(manifest_ch)

    UCSC_BLAT(MANIFEST2FASTA.out.fasta, genome_ch)

    ch_versions = ch_versions.mix(UCSC_BLAT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
