#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:56:04 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import re
import logging
import argparse

import Bio.SeqIO
import Bio.SearchIO

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Parse chip alignment file")
parser.add_argument(
    "-a", "--alignment", required=True, help="Chip alignment file")
parser.add_argument(
    "-c", "--chip_sequence", required=True, help="Chip fasta file")
parser.add_argument(
    "-g", "--genome_sequence", required=True, help="Genome fasta file")
args = parser.parse_args()

pattern = re.compile(
    ".* snp (.*), iln_snp: (.*), iln_pos: (.*), iln_strand: (.*), "
    "src_strand: (.*)")


def parse_description(description):
    match = re.search(pattern, description)
    snp, iln_snp, iln_pos, iln_strand, src_strand = match.groups()
    return snp, iln_snp, int(iln_pos), iln_strand, src_strand


def check_strand(hsp):
    if hsp.hit_strand < 0:
        raise NotImplementedError("Reversed hit strand not yet managed")

    if hsp.query_strand > 0:
        return "forward"

    else:
        return "reverse"


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)

    # read the chip fasta file
    chip_sequences = Bio.SeqIO.index(args.chip_sequence, "fasta")

    # read the genome file
    genome_sequences = Bio.SeqIO.index(args.genome_sequence, "fasta")

    for result in Bio.SearchIO.parse(args.alignment, "blat-psl", pslx=True):
        logger.debug("-------------------------------------------------------")
        # result represent a single search query

        # collect SNP info
        snp_sequence = chip_sequences[result.id]
        snp, iln_snp, iln_pos, iln_strand, src_strand = parse_description(
            snp_sequence.description)

        logger.debug(
            f"Processing SNP: {result.id}, snp: {snp}, iln_snp: {iln_snp}, "
            f"iln_pos: {iln_pos}, iln_strand: {iln_strand}, "
            f"src_strand: {src_strand}")

        # hit represents a single database hit
        logger.debug(f"Got {len(result.hits)} hits for {result.id}")

        for hit in result.hits:
            logger.debug(f"Processing hit {hit.id} for {result.id}")

            chr_sequence = genome_sequences[hit.id]

            # hsp represents region(s) of significant alignments between
            # query and hit sequences
            logger.debug(f"Got {len(hit.hsps)} hsp for {hit.id}")

            for hsp in hit.hsps:
                logger.debug(
                    f"Query range {hsp.query_range} ({hsp.query_strand}) "
                    f"Hit range {hsp.hit_range} ({hsp.hit_strand})")

                orient = check_strand(hsp)

                logger.debug(f"Orient is '{orient}'")

                # get REF sequence from reference genome
                ref_pos = hsp.hit_start + iln_pos

                # this is 0-based index
                ref_allele = chr_sequence[ref_pos-1]

                logger.debug(
                    f"Reference allele: {ref_allele} at "
                    f"{hit.id}:{ref_pos}"
                )

                # TODO: define ALT allele

        # debug
        # break
