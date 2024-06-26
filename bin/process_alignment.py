#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 16:56:04 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import csv
import logging
import argparse

import Bio.SeqIO
import Bio.SearchIO
import Bio.AlignIO

from helper.blast import BlastResult

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Parse chip alignment file",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-a", "--alignment", required=True, help="Chip alignment file")
parser.add_argument(
    "--output_aln", required=False, help="Output alignment file")
parser.add_argument(
    "-o", "--output_csv", required=True, help="Output CSV file")
parser.add_argument(
    "--error_csv", required=False,
    help="SNP which can't be mapped for any reason")
parser.add_argument(
    "--length_pct", required=False, type=float, default=95,
    help="Percentage of the query aligned (default: %(default)s)")
parser.add_argument(
    "--ident_pct", required=False, type=float, default=97,
    help="Percentage identities in alignment (default: %(default)s)")


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # process received arguments
    args = parser.parse_args()

    # open file for writing
    output_csv_fh = open(args.output_csv, "w")
    writer = csv.writer(output_csv_fh, delimiter=",", lineterminator="\n")

    header = [
        "snp_name", "chrom", "position", "alleles", "illumina",
        "illumina_forward", "illumina_strand", "strand",
        "ref", "alt"
    ]

    writer.writerow(header)

    if args.output_aln:
        output_aln_fh = open(args.output_aln, "w")

    # track errors in a file
    if args.error_csv:
        error_csv_fh = open(args.error_csv, "w")
        error = csv.writer(error_csv_fh, delimiter=",", lineterminator="\n")
        error.writerow(["snp_name", "illumina", "illumina_strand", "reason"])

    logger.info(f"Processing {args.alignment}")
    logger.info(
        f"Filtering for length_pct: {args.length_pct}% and ident_pct: "
        f"{args.ident_pct}%")

    for record in Bio.SearchIO.parse(args.alignment, "blast-xml"):
        # result represent a single search query
        logger.info("-------------------------------------------------------")
        result = BlastResult(record)

        # filter alignments
        result.filter_results(
            length_pct=args.length_pct, ident_pct=args.ident_pct)

        # process SNP informations
        lines, alignments, discarded_snps = result.process_alignments()

        for line in lines:
            logger.debug(f"Writing {line}")
            writer.writerow(line)

        # write alignments
        if args.output_aln:
            for aln in alignments:
                Bio.AlignIO.write([aln], output_aln_fh, "clustal")

        # write errors
        if args.error_csv:
            for discarded in discarded_snps:
                error.writerow(discarded)

        # debug
        # break

    output_csv_fh.close()

    if args.output_aln:
        output_aln_fh.close()

    if args.error_csv:
        error_csv_fh.close()

    logger.info("-------------------------------------------------------")
    logger.info("Done!")
