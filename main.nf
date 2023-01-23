#!/usr/bin/env nextflow

// DSL 2
nextflow.enable.dsl = 2

// include workflow dependencies from external modules
include { MANIFEST2FASTA } from './modules/local/manifest2fasta'
include { TABIX_BGZIP } from './modules/nf-core/tabix/bgzip/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN } from './modules/nf-core/blast/blastn/main'
include { PROCESSALIGNMENT } from './modules/local/processalignment'
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

    MANIFEST2FASTA.out.fasta.map{
        fasta -> [[ id: fasta.baseName ], fasta]
    }.set{ blast_input }

    BLAST_BLASTN(blast_input, BLAST_MAKEBLASTDB.out.db)
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    PROCESSALIGNMENT(BLAST_BLASTN.out.txt)

    PROCESSALIGNMENT.out.csv.collectFile(
        name: "results.csv",
        storeDir: "${params.outdir}/chip_aligned",
        keepHeader: true,
    )

    PROCESSALIGNMENT.out.aln.collectFile(
        name: "results.aln",
        storeDir: "${params.outdir}/chip_aligned",
        keepHeader: false,
    )

    PROCESSALIGNMENT.out.err.collectFile(
        name: "results.err",
        storeDir: "${params.outdir}/chip_aligned",
        keepHeader: true,
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
