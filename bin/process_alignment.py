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

from helper.blat import BlatResult
from helper.utils import text_or_gzip_open

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Parse chip alignment file",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
parser.add_argument(
    "--lenth_pct", required=False, type=float, default=95,
    help="Percentage of the query aligned")
parser.add_argument(
    "--ident_pct", required=False, type=float, default=97,
    help="Percentage identities in alignment")
args = parser.parse_args()


def parse_chromosome(sequence):
    # try to determin chromosome name
    match_chrom = re.search("chromosome (.*),", sequence.description)
    match_scaff = re.search("(scaffold_[0-9]+),", sequence.description)
    match_unknw = re.search("(unplaced_[0-9]+),", sequence.description)

    chrom = sequence.id

    if match_chrom:
        chrom = match_chrom.groups()[0]

    elif match_scaff:
        chrom = match_scaff.groups()[0]

    elif match_unknw:
        chrom = match_unknw.groups()[0]

    return chrom


def parse_chromosomes(genome_file):
    id2chromosome = {}

    with text_or_gzip_open(genome_file) as handle:
        for sequence in Bio.SeqIO.parse(handle, "fasta"):
            id2chromosome[sequence.id] = parse_chromosome(sequence)

    return id2chromosome


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # read the chip fasta file
    chip_sequences = Bio.SeqIO.index(args.chip_sequence, "fasta")

    # get descriptions from genome file
    id2chromosome = parse_chromosomes(args.genome_sequence)

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
        # result represent a single search query
        logger.info("-------------------------------------------------------")
        result = BlatResult(result)

        # collect SNP info
        snp_sequence = chip_sequences[result.snp_id]
        result.read_sequence_manifest(snp_sequence)

        # filter alignments
        result.filter_results(
            lenth_pct=args.lenth_pct, ident_pct=args.ident_pct)

        # process SNP informations
        lines, alignments = result.process_alignments(id2chromosome)

        for line in lines:
            logger.debug(f"Writing {line}")
            writer.writerow(line)

        # write alignments
        if args.output_aln:
            for aln in alignments:
                Bio.AlignIO.write([aln], output_aln_fh, "clustal")

        # debug
        # break

    output_csv_fh.close()

    if args.output_aln:
        output_aln_fh.close()
