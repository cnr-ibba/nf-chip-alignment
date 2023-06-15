#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { MANIFEST2FASTA } from './modules/local/manifest2fasta'
include { TABIX_BGZIP } from './modules/nf-core/tabix/bgzip/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_MEGABLAST } from './modules/nf-core/blast/blastn/main'
include { PROCESSALIGNMENT } from './modules/local/processalignment'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'


workflow {
    // get manifest from parameters
    manifest_ch = Channel.fromPath( params.manifest )
    genome_ch = Channel.fromPath( params.genome )
        .map{ it -> [[id:it.baseName], it] }
    ch_versions = Channel.empty()

    if (params.genome.endsWith('.gz')) {
        // unpack genome
        TABIX_BGZIP(genome_ch)
        genome_ch = TABIX_BGZIP.out.output
    }

    BLAST_MAKEBLASTDB(genome_ch)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // convert a manifest into fasta
    MANIFEST2FASTA(manifest_ch)

    MANIFEST2FASTA.out.fasta.map{
        fasta -> [[ id: fasta.baseName ], fasta]
    }.set{ blast_input }

    BLAST_MEGABLAST(blast_input, BLAST_MAKEBLASTDB.out.db)
    ch_versions = ch_versions.mix(BLAST_MEGABLAST.out.versions)

    BLAST_MEGABLAST.out.txt.map{
        meta, alignment -> [[ id: alignment.baseName ], alignment]
    }.set{ processalignment_input }

    PROCESSALIGNMENT(processalignment_input)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
