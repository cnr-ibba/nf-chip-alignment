#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { MANIFEST2FASTA } from './modules/local/manifest2fasta'
include { UCSC_FATOTWOBIT } from './modules/local/ucsc/fatotwobit'
include { UCSC_BLAT } from './modules/local/ucsc/blat'
include { PROCESSALIGNMENT } from './modules/local/processalignment'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/modules/custom/dumpsoftwareversions/main'

// get manifest from parameters
manifest_ch = Channel.fromPath( params.manifest )
genome_ch = Channel.fromPath( params.genome )
ch_versions = Channel.empty()


workflow {
    // convert a manifest into fasta
    MANIFEST2FASTA(manifest_ch)

    // database preparation
    UCSC_FATOTWOBIT(genome_ch)

    // probe alignment
    UCSC_BLAT(MANIFEST2FASTA.out.fasta, UCSC_FATOTWOBIT.out.twobit)

    PROCESSALIGNMENT(MANIFEST2FASTA.out.fasta, genome_ch, UCSC_BLAT.out.pslx)

    ch_versions = ch_versions.mix(UCSC_BLAT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
