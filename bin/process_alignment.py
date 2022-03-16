#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:56:04 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import re
import csv
import logging
import argparse

import Bio.SeqIO
import Bio.SearchIO
import Bio.AlignIO

from helper import complement

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Parse chip alignment file")
parser.add_argument(
    "-a", "--alignment", required=True, help="Chip alignment file")
parser.add_argument(
    "-c", "--chip_sequence", required=True, help="Chip fasta file")
parser.add_argument(
    "-g", "--genome_sequence", required=True, help="Genome fasta file")
parser.add_argument(
    "--output_aln", required=False, help="Output alignment file")
parser.add_argument(
    "-o", "--output_csv", required=True, help="Output CSV file")
args = parser.parse_args()

pattern = re.compile(
    ".* iln_snp: (.*), iln_pos: (.*), iln_strand: (.*)")


def parse_description(description):
    match = re.search(pattern, description)
    iln_snp, iln_pos, iln_strand = match.groups()
    return iln_snp, int(iln_pos), iln_strand


def check_strand(hsp):
    if hsp.hit_strand < 0:
        raise NotImplementedError("Reversed hit strand not yet managed")

    if hsp.query_strand > 0:
        return "forward"

    else:
        return "reverse"


def get_alt_allele(snp, ref, orient):
    """Get ALT allele giving the ref and orientation"""

    if orient == 'reverse':
        snp = complement(snp)

    genotype = snp.split("/")

    ref_index = genotype.index(ref)

    # the alt allele is the non reference allele
    if ref_index == 0:
        return genotype[1]
    else:
        return genotype[0]


def get_illumina_forward(snp, orient):
    """Get the illumina forward SNP. If the strand is forward, is the same of
    received from illumina. If not, reverse the SNP"""

    if orient == 'reverse':
        snp = complement(snp)

    return snp


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)

    # read the chip fasta file
    chip_sequences = Bio.SeqIO.index(args.chip_sequence, "fasta")

    # read the genome file
    genome_sequences = Bio.SeqIO.index(args.genome_sequence, "fasta")

    # open file for writing
    output_csv_fh = open(args.output_csv, "w")
    writer = csv.writer(output_csv_fh, delimiter=",", lineterminator="\n")

    if args.output_aln:
        output_aln_fh = open(args.output_aln, "w")

    header = [
        "snp_name", "chrom", "position", "alleles", "illumina",
        "illumina_forward", "illumina_strand", "strand",
        "ref", "alt"
    ]

    writer.writerow(header)

    for result in Bio.SearchIO.parse(args.alignment, "blat-psl", pslx=True):
        logger.debug("-------------------------------------------------------")
        # result represent a single search query

        # collect SNP info
        snp_sequence = chip_sequences[result.id]
        iln_snp, iln_pos, iln_strand = parse_description(
            snp_sequence.description)

        logger.debug(
            f"Processing SNP: {result.id}, iln_snp: {iln_snp}, "
            f"iln_pos: {iln_pos}, iln_strand: {iln_strand}")

        # hit represents a single database hit
        logger.debug(f"Got {len(result.hits)} hits for {result.id}")

        if len(result.hits) > 1:
            raise NotImplementedError("More than one hits!")

        for hit in result.hits:
            logger.debug(f"Processing hit {hit.id} for {result.id}")

            chr_sequence = genome_sequences[hit.id]

            # hsp represents region(s) of significant alignments between
            # query and hit sequences
            logger.debug(f"Got {len(hit.hsps)} hsp for {hit.id}")

            if len(hit.hsps) > 1:
                raise NotImplementedError("More than one hsp!")

            for hsp in hit.hsps:
                logger.debug(
                    f"Query range {hsp.query_range} ({hsp.query_strand}), "
                    f"Hit range {hsp.hit_range} ({hsp.hit_strand}), "
                    f"Score {hsp.score}, ident_pct {hsp.ident_pct}")

                orient = check_strand(hsp)

                logger.debug(f"Orient is '{orient}'")

                # get REF sequence from reference genome

                # find 'n' position on query (where SNP is in alignment)
                # is a 0 based position
                aln_pos = hsp.query.seq.find('n')

                ref_pos = hsp.hit_start + aln_pos

                # this is 0-based index
                ref_allele = chr_sequence[ref_pos]

                logger.debug(
                    f"Reference allele: {ref_allele} at "
                    f"{hit.id}:{ref_pos+1}"
                )

                alt_allele = get_alt_allele(iln_snp, ref_allele, orient)

                logger.debug(
                    f"Alternative allele: {alt_allele} at "
                    f"{hit.id}:{ref_pos+1}"
                )

                # alleles a la NCBI
                alleles = "/".join(sorted([ref_allele, alt_allele]))

                # define a new data row. Mind to 1-based positions
                line = [
                    result.id, hit.id, ref_pos+1, alleles, iln_snp,
                    get_illumina_forward(iln_snp, orient), iln_strand, orient,
                    ref_allele, alt_allele]
                writer.writerow(line)

                # annotate alignment
                hsp.aln[-1].id = "{0}:{1}-{2}".format(
                    hsp.aln[1].id.split("|")[-1],
                    hsp.hit_range[0],
                    hsp.hit_range[1])

                # write alignment
                if args.output_aln:
                    Bio.AlignIO.write([hsp.aln], output_aln_fh, "clustal")

        # debug
        # break

    output_csv_fh.close()

    if args.output_aln:
        output_aln_fh.close()
