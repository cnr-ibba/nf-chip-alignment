#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:59:36 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import re
import logging

from Bio.SearchIO._model.query import QueryResult
from Bio.SeqRecord import SeqRecord

from .utils import complement

# Get an instance of a logger
logger = logging.getLogger(__name__)

manifest_pattern = re.compile(
    ".* iln_snp: (.*), iln_pos: (.*), iln_strand: (.*)")


def parse_description(description):
    match = re.search(manifest_pattern, description)
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


class BlatResult():
    queryresult = None
    filtered = None
    is_filtered = False
    status = None
    snp_id = None
    iln_snp = None
    iln_pos = None
    iln_strand = None
    probe_len = None
    best_hit = None

    def __init__(self, queryresult: QueryResult = None):
        if queryresult:
            self.queryresult = queryresult
            self.snp_id = queryresult.id

    def __repr__(self):
        return (
            f"BlatResult: {self.snp_id}, iln_snp: {self.iln_snp}, "
            f"iln_pos: {self.iln_pos}, iln_strand: {self.iln_strand} "
            f"probe length {self.probe_len}")

    def read_sequence_manifest(self, sequence: SeqRecord):
        # collect SNP info
        self.iln_snp, self.iln_pos, self.iln_strand = parse_description(
            sequence.description)
        self.probe_len = self.queryresult.seq_len

        logger.info(
            f"Processing SNP: {self.snp_id}, iln_snp: {self.iln_snp}, "
            f"iln_pos: {self.iln_pos}, iln_strand: {self.iln_strand} "
            f"probe length {self.probe_len}")

    def filter_results(self, lenth_pct=95, ident_pct=97):
        # reset best hsp
        self.best_hit = None

        # hit represents a single database hit
        logger.debug(
            f"Got {len(self.queryresult.hits)} hits for {self.queryresult.id}")

        # filter results by score (query aligned)
        def filter_hsps(hsp):
            if hsp.is_fragmented:
                logger.debug(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Found {len(hsp.fragments)} fragments"
                )
                return False

            # score is a function of probe length
            if (hsp.score < self.probe_len * lenth_pct / 100 or
                    hsp.ident_pct < ident_pct):
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

        filtered = self.queryresult.hsp_filter(filter_hsps)

        logger.info(f"Got {len(filtered.hits)} hits after filtering")

        if len(filtered.hits) == 0 or len(filtered.hsps) == 0:
            logger.warning(
                f"All alignments have been filtered out for "
                f"{self.queryresult.id}")

            self.status = "No valid alignments after filtering"

        elif len(filtered.hits) > 1 or len(filtered.hsps) > 1:
            for hsp in filtered.hsps:
                logger.debug(
                    f"{hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Score: {hsp.score} (ident_pct {hsp.ident_pct})"
                )
            logger.warning(
                f"Got {len(filtered.hsps)} alignments after filtering for "
                f"{self.queryresult.id}")

            self.status = "Too many alignments after filtering"

        else:
            self.best_hit = filtered[0]
            self.status = "Found valid alignment after filtering"

        self.filtered = filtered
        self.is_filtered = True

    def process_alignments(self, id2chromosome):
        """Returns an output record and the processed alignment"""

        self.lines, self.alignments, discarded_snps = [], [], []

        if self.is_filtered and not self.best_hit:
            logger.warning(f"Discarding {self}")
            # return an empty row
            line = [
                self.snp_id, 0, 0, None, self.iln_snp, None,
                self.iln_strand, None, None, None]
            self.lines = [line]

            discarded_snps = [[
                self.snp_id,
                self.iln_snp,
                self.iln_strand,
                self.status]]

        else:
            for line, alignment, discarded in self.__process_hits(
                    id2chromosome):
                self.lines.append(line)
                self.alignments.append(alignment)

                if discarded:
                    discarded_snps.append(discarded)

        return self.lines, self.alignments, discarded_snps

    def __process_hits(self, id2chromosome):
        discarded = []

        for i, hit in enumerate(self.filtered.hits):
            logger.info(f"Processing hit {i}: {hit.id} for {self.snp_id}")

            # attempt to determine chromosome name
            chrom = id2chromosome[hit.id]

            logger.info(f"Detected chromosome for {hit.id} is '{chrom}'")

            # hsp represents region(s) of significant alignments between
            # query and hit sequences
            logger.debug(f"Got {len(hit.hsps)} hsp for {hit.id}")

            for j, hsp in enumerate(hit.hsps):
                logger.info(
                    f"Hsp {j}: has {len(hsp.fragments)} fragments. "
                    f"Query range {hsp.query_range_all} "
                    f"({hsp.query_strand_all}), "
                    f"Hit range {hsp.hit_range_all} ({hsp.hit_strand_all}), "
                    f"Score {hsp.score}, ident_pct {hsp.ident_pct}")

                orient = check_strand(hsp)

                logger.info(f"Orient is '{orient}'")

                # the SNP position in the alignment, supposing no gap in
                # query sequence (mind to the query strand)
                if hsp.query_strand > 0:
                    aln_pos = self.iln_pos - hsp.query_start - 1
                else:
                    aln_pos = hsp.query_end - self.iln_pos

                # get and annotate alignment
                alignment = hsp.aln
                alignment[-1].id = "{0}:{1}-{2}".format(
                    chrom,
                    hsp.hit_range[0],
                    hsp.hit_range[1])

                # check that is letter is a N
                if alignment[0][aln_pos].upper() != 'N':
                    logger.error(alignment[:, aln_pos-5:aln_pos+6])
                    raise Exception(
                        f"Cannot find the SNP in position {aln_pos}")

                ref_pos = hsp.hit_start + aln_pos

                # this is 0-based index
                ref_allele = alignment[1][aln_pos].upper()

                logger.info(
                    f"Reference allele: {ref_allele} at "
                    f"{hit.id}:{ref_pos+1}"
                )

                try:
                    alt_allele = get_alt_allele(
                        self.iln_snp, ref_allele, orient)

                except ValueError:
                    logger.warning(
                        f"Cannot find alt_allele in illumina SNP: "
                        f"'{ref_allele}' not in {self.iln_snp}")
                    logger.warning(f"Discarding {self}")
                    self.status = "Allele doesn't match to reference"

                    line = [
                        self.snp_id, 0, 0, None, self.iln_snp, None,
                        self.iln_strand, None, None, None]
                    discarded = [
                        self.snp_id,
                        self.iln_snp,
                        self.iln_strand,
                        self.status]

                else:
                    logger.info(
                        f"Alternative allele: {alt_allele} at "
                        f"{hit.id}:{ref_pos+1}"
                    )

                    # alleles a la NCBI
                    alleles = "/".join(sorted([ref_allele, alt_allele]))

                    # define a new data row. Mind to 1-based positions
                    line = [
                        self.snp_id, chrom, ref_pos+1, alleles, self.iln_snp,
                        get_illumina_forward(self.iln_snp, orient),
                        self.iln_strand, orient, ref_allele, alt_allele]

                yield line, alignment, discarded
