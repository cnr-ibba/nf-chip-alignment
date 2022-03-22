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

        seq = Bio.Seq.Seq(record.sourceseq)

        # read snp as illuSNP
        try:
            iln_snp = IlluSNP(record.sourceseq, max_iter=60)

        except IlluSNPException as exc:
            logger.warning(exc)
            logger.warning(
                f"Can't determine illumina strand for '{record.name}'")
            continue

        logger.debug(
            f"id: '{record.name}', iln_snp: {iln_snp.snp}, iln_strand: "
            f"{iln_snp.strand}")

        # logger.debug(f"topgenomicseq: {record.topgenomicseq}")
        logger.debug(f"record.sourceseq: {seq}")

        # find SNP in sequence
        start, end = iln_snp.pos

        logger.debug(f"found {iln_snp.snp} in {start}:{end}")

        # trasnform in a mutable object
        seq = seq.tomutable()

        # replace SNP with ambiguos code
        seq[start:end] = SNP2BASES[iln_snp.snp]

        record = Bio.SeqRecord.SeqRecord(
            id=record.name,
            name=record.name,
            description=(
                f"iln_snp: {iln_snp.snp}, iln_pos: {start+1}, "
                f"iln_strand: {iln_snp.strand}"
            ),
            seq=seq)

        Bio.SeqIO.write([record], handle, "fasta")

    handle.close()
    logger.info("Done!")
