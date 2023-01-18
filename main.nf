#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { MANIFEST2FASTA } from './modules/local/manifest2fasta'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN } from './modules/nf-core/blast/blastn/main'
// include { PROCESSALIGNMENT } from './modules/local/processalignment'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'

// get manifest from parameters
manifest_ch = Channel.fromPath( params.manifest )
genome_ch = Channel.fromPath( params.genome )
ch_versions = Channel.empty()


workflow {
    // convert a manifest into fasta
    MANIFEST2FASTA(manifest_ch)

    BLAST_MAKEBLASTDB(genome_ch)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    MANIFEST2FASTA.out.fasta.map{
        fasta -> [[ id: fasta.baseName ], fasta]
    }.set{ blast_input }

    BLAST_BLASTN(blast_input, BLAST_MAKEBLASTDB.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    // PROCESSALIGNMENT(MANIFEST2FASTA.out.fasta, genome_ch, UCSC_BLAT.out.pslx)

    // ch_versions = ch_versions.mix(UCSC_BLAT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
