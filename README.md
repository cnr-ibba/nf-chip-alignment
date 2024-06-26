
# Nextflow chip alignment

This nextflow pipeline is an attempt to align probes from a chip-manifest file
on a reference genomes using [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
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

## About alignments

Alignments are made using `megablast` task, which is more stringent and faster
then the default `blastn`.
Alignments are filtered relying *percentage aligned* and *percentage
identity*. A probe is considered unmapped if all alignments are filtered,
if there are more than one alignment after filtering or if the SNP doesn't
match to the reference genome. The last condition however is not correlated in
a issue in alignment, however it could be difficult to join this particular SNP
with other data or refer it to a reference SNP. For such reason SNP is discarded
even if there's a perfect match of the probe
