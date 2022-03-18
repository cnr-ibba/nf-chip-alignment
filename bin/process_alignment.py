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


def parse_chromosome(sequence):
    # try to determin chromosome name
    match_chrom = re.search("chromosome (.*),", sequence.description)
    match_scaff = re.search("(scaffold_[0-9]+),", sequence.description)
    match_unknw = re.search("(unplaced_[0-9]+),", sequence.description)

    chrom = hit.id

    if match_chrom:
        chrom = match_chrom.groups()[0]

    elif match_scaff:
        chrom = match_scaff.groups()[0]

    elif match_unknw:
        chrom = match_unknw.groups()[0]

    return chrom


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
            f"iln_pos: {iln_pos}, iln_strand: {iln_strand} "
            f"probe length {result.seq_len}")

        # hit represents a single database hit
        logger.debug(f"Got {len(result.hits)} hits for {result.id}")

        # filter results by score (query aligned)
        def filter_hsps(hsp):
            if hsp.is_fragmented:
                logger.debug(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Found {len(hsp.fragments)} fragments"
                )
                return False

            if hsp.score < result.seq_len * 0.95 or hsp.ident_pct < 97:
                logger.debug(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Bad Score: {hsp.score} (ident_pct {hsp.ident_pct})"
                )
                return False

            logger.debug(
                f"Keeping {hsp.hit_id}:{hsp.hit_range_all}: "
                f"Score: {hsp.score} (ident_pct {hsp.ident_pct})"
            )

            return True

        filtered = result.hsp_filter(filter_hsps)

        logger.debug(f"Got {len(filtered.hits)} hits after filtering")

        if len(filtered.hits) == 0 or len(filtered.hsps) == 0:
            logger.warning(
                f"All alignements have been filtered out for {result.id}")

            # write an empty row
            line = [
                result.id, 0, 0, None, iln_snp, None, iln_strand, None,
                None, None]
            logger.debug(f"Writing {line}")
            writer.writerow(line)

            # skip processing
            continue

        elif len(filtered.hits) > 1 or len(filtered.hsps) > 1:
            for hsp in filtered.hsps:
                logger.warning(
                    f"{hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Score: {hsp.score} (ident_pct {hsp.ident_pct})"
                )
            logger.error("Got more alignment than expected for {result.id}")

            # write an empty row
            line = [
                result.id, 0, 0, None, iln_snp, None, iln_strand, None,
                None, None]
            logger.debug(f"Writing {line}")
            writer.writerow(line)

            # skip processing
            continue

        for i, hit in enumerate(filtered.hits):
            logger.debug(f"Processing hit {i}: {hit.id} for {result.id}")

            chr_sequence = genome_sequences[hit.id]

            # attempt to determine chromosome name
            chrom = parse_chromosome(chr_sequence)

            logger.debug(f"Detected chromosome for {hit.id} is {chrom}")

            # hsp represents region(s) of significant alignments between
            # query and hit sequences
            logger.debug(f"Got {len(hit.hsps)} hsp for {hit.id}")

            for j, hsp in enumerate(hit.hsps):
                logger.debug(
                    f"Hsp {j}: has {len(hsp.fragments)} fragments. "
                    f"Query range {hsp.query_range_all} "
                    f"({hsp.query_strand_all}), "
                    f"Hit range {hsp.hit_range_all} ({hsp.hit_strand_all}), "
                    f"Score {hsp.score}, ident_pct {hsp.ident_pct}")

                orient = check_strand(hsp)

                logger.debug(f"Orient is '{orient}'")

                # get REF sequence from reference genome

                # find 'n' position on query (where SNP is in alignment)
                # is a 0 based position
                aln_pos = hsp.query.seq.find('n')

                ref_pos = hsp.hit_start + aln_pos

                # this is 0-based index
                ref_allele = chr_sequence[ref_pos].upper()

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
                    result.id, chrom, ref_pos+1, alleles, iln_snp,
                    get_illumina_forward(iln_snp, orient), iln_strand, orient,
                    ref_allele, alt_allele]
                logger.debug(f"Writing {line}")
                writer.writerow(line)

                # annotate alignment
                hsp.aln[-1].id = "{0}:{1}-{2}".format(
                    chrom,
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
