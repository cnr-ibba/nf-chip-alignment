#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 12:33:06 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import logging
import argparse

import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord

from helper import SNP2BASES
from helper.illumina import read_Manifest, IlluSNP, IlluSNPException

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Extract fasta sequences from chip manifest file")
parser.add_argument(
    "-i", "--input", required=True, help="Manifest file input path")
parser.add_argument(
    "-o", "--output", required=True, help="Fasta file output path")
args = parser.parse_args()


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    logger.info(f"Opening '{args.output}' for writing")
    handle = open(args.output, "w")

    for i, record in enumerate(read_Manifest(args.input, delimiter=",")):
        logger.debug("-------------------------------------------------------")

        # create a new sequence object by strippint first and last Ns
        sequence = Bio.Seq.Seq(record.sourceseq).strip(chars="N")
        logger.debug(f"Got sequence: {sequence}")

        # read snp as illuSNP
        try:
            iln_snp = IlluSNP(str(sequence), max_iter=60)

        except IlluSNPException as exc:
            logger.warning(exc)
            logger.warning(
                f"Can't determine illumina strand for '{record.name}'")
            continue

        logger.debug(
            f"id: '{record.name}', iln_snp: {iln_snp.illumina}, iln_strand: "
            f"{iln_snp.strand}")

        # find SNP in sequence
        start, end = iln_snp.pos

        # this are the alleles in sequence fomatted like NCBI
        logger.debug(f"found {iln_snp.alleles} in {start}:{end}")

        # transform in a mutable object
        sequence = sequence.tomutable()

        # replace SNP with ambiguous code
        sequence[start:end] = SNP2BASES[iln_snp.alleles]

        record = Bio.SeqRecord.SeqRecord(
            id=record.name,
            name=record.name,
            description=(
                f"iln_snp: {iln_snp.illumina}, iln_pos: {start+1}, "
                f"iln_strand: {iln_snp.strand}"
            ),
            seq=sequence)

        Bio.SeqIO.write([record], handle, "fasta")

    handle.close()
    logger.info("Done!")
