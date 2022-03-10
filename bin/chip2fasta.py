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
from helper.illumina import read_Manifest

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
    logging.basicConfig(level=logging.DEBUG, format=log_fmt)

    logger.info(f"Opening '{args.output}' for writing")
    handle = open(args.output, "w")

    for i, record in enumerate(read_Manifest(args.input, delimiter=",")):
        seq = Bio.Seq.MutableSeq(record.sourceseq)

        logger.debug(
            f"iln_snp: {record.snp}, iln_strand: {record.ilmnstrand} "
            f"src_strand: {record.sourcestrand}")
        logger.debug(
            f"{seq}")

        # find SNP in sequence
        idx = seq.find(record.snp)

        if idx < 0:
            raise Exception(f"Snp '[{record.snp}]' not found in {seq}")

        # get snp interval. Mind to []
        start = idx - 1
        end = idx + len(record.snp) + 1

        logger.debug(f"found {record.snp} in {start}:{end}")

        # replace SNP with ambiguos code
        seq[start:end] = SNP2BASES[record.snp]

        record = Bio.SeqRecord.SeqRecord(
            id=record.name,
            name=record.name,
            description=(
                f"iln_snp: {record.snp}, iln_strand: {record.ilmnstrand} "
                f"src_strand: {record.sourcestrand}"
            ),
            seq=seq)

        Bio.SeqIO.write([record], handle, "fasta")

    handle.close()
    logger.info("Done!")
