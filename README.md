
# Nextflow chip alignment

This nextflow pipeline is an attempt to align probes from a chip-manifest file
on a reference genomes using [blat](https://genome.ucsc.edu/goldenPath/help/blatSpec.html).
The reason of doing this is to define the strand orientation of the probe on the
reference sequence, in order to be able to convert a genotype from *illumina top*
to *forward* (or *ALT/REF*) coordinates. This pipeline was developed as a component
of the [SMARTER project](https://www.smarterproject.eu/) (WP4 in particular) and
will be the first step towards a new [SNPchiMp](https://webserver.ibba.cnr.it/SNPchimp/)
release.

## Requisites

To run this nextflow pipeline, you require both a genome reference sequence (better
to download it from NCBI ftp) and a manifest *csv* file from an *illumina/affymetrix*
chip. Next, you you will need `nextflow` installed and at least one of this
different executors: `conda`, `singularity` and `docker`. You can choose to clone
this repository if you plan to change this pipeline according your needs:

```bash
git clone https://github.com/cnr-ibba/nf-chip-alignment.git
```

The other way to running this pipeline is described in
[pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html#pipeline-sharing)
nextflow manual.

## Testing pipeline

You can test your local configuration by calling the `test` profile with one of
the three supported executors. For example, to test this pipeline with `singularity`:

```bash
nextflow run cnr-ibba/nf-chip-alignment -resume -profile test,singularity
```

## Running pipeline

You can call this pipeline with you local data using the `--genome` and `--manifest`
options, for example:

```bash
nextflow run cnr-ibba/nf-chip-alignment -resume -profile singularity --manifest <manifest.csv> --genome <reference genome>
```

*gzipped* input files are supported.
