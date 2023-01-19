#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { MANIFEST2FASTA } from './modules/local/manifest2fasta'
include { TABIX_BGZIP } from './modules/nf-core/tabix/bgzip/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN } from './modules/nf-core/blast/blastn/main'
// include { PROCESSALIGNMENT } from './modules/local/processalignment'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'

// get manifest from parameters
manifest_ch = Channel.fromPath( params.manifest )
genome_ch = Channel.fromPath( params.genome )
ch_versions = Channel.empty()


workflow {
    if (params.genome.endsWith('.gz')) {
        // unpack genome
        TABIX_BGZIP(genome_ch.map{ it -> [[id:it.baseName], it] })
        genome_ch = TABIX_BGZIP.out.output.map{ meta, fasta -> [fasta] }
    }

    BLAST_MAKEBLASTDB(genome_ch)
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    // convert a manifest into fasta
    MANIFEST2FASTA(manifest_ch)

    // split fasta input file in chunks
    MANIFEST2FASTA.out.fasta
        .splitFasta( by: 1000, file: true, elem: 1 )
        .map {
            fasta -> [[ id: fasta.baseName ], fasta ]
        }
        .set { chunks }

    chunks.combine(BLAST_MAKEBLASTDB.out.db)
        .multiMap { record ->
            chunk: [record[0], record[1]]
            blast_db: record[2]
        }.set {
            blast_input
        }

    BLAST_BLASTN(blast_input.chunk, blast_input.blast_db)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    // PROCESSALIGNMENT(MANIFEST2FASTA.out.fasta, genome_ch, UCSC_BLAT.out.pslx)

    // ch_versions = ch_versions.mix(UCSC_BLAT.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
