#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 13:59:36 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import re
import logging
import Bio.SeqIO

from Bio.SearchIO._model.hsp import HSP
from Bio.SearchIO._model.query import QueryResult
from Bio.SeqRecord import SeqRecord

from .utils import complement, text_or_gzip_open

# Get an instance of a logger
logger = logging.getLogger(__name__)

manifest_pattern = re.compile(
    ".* iln_snp: (.*), iln_pos: (.*), iln_strand: (.*)")


def parse_chromosomes(genome_file):
    id2chromosome = {}

    def parse_chromosome(sequence):
        # try to determine chromosome name
        match_chrom = re.search(r"chromosome (\w*),", sequence.description)
        match_scaff = re.search(r"(scaffold_[0-9]+),", sequence.description)
        match_unknw = re.search(r"(unplaced_[0-9]+),", sequence.description)

        chrom = sequence.id

        if match_chrom:
            chrom = match_chrom.groups()[0]

        elif match_scaff:
            chrom = match_scaff.groups()[0]

        elif match_unknw:
            chrom = match_unknw.groups()[0]

        return chrom

    with text_or_gzip_open(genome_file) as handle:
        for sequence in Bio.SeqIO.parse(handle, "fasta"):
            id2chromosome[sequence.id] = parse_chromosome(sequence)

    return id2chromosome


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


def get_alt_allele(snp, ref, hsp):
    """Get ALT allele giving the ref and HSP"""

    orient = check_strand(hsp)

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


def annotate_alignment(hsp: HSP, chrom: str):
    """Annotate alignment with positions and chromosome names"""

    alignment = hsp.aln
    alignment[0].id = "{0}:{1}-{2}".format(
        alignment[0].id,
        hsp.query_range[0],
        hsp.query_range[1])
    alignment[1].id = "{0}:{1}-{2}".format(
        chrom,
        hsp.hit_range[0],
        hsp.hit_range[1])

    return alignment


class BlatException(Exception):
    """Base exception class for BlatResult"""

    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        return repr(self.value)


class BlatResult():
    queryresult = None
    filtered = None
    is_filtered = False
    got_manifest = False
    status = None
    snp_id = None
    iln_snp = None
    iln_pos = None  # 1 based position
    iln_strand = None
    probe_len = None
    best_hit = None

    def __init__(self, queryresult: QueryResult = None):
        if queryresult:
            self.queryresult = queryresult
            self.snp_id = queryresult.id
            self.probe_len = self.queryresult.seq_len

    def __repr__(self):
        return (
            f"BlatResult: {self.snp_id}, iln_snp: {self.iln_snp}, "
            f"iln_pos: {self.iln_pos}, iln_strand: {self.iln_strand}, "
            f"probe_len: {self.probe_len}")

    def read_sequence_manifest(self, sequence: SeqRecord):
        # collect SNP info
        self.iln_snp, self.iln_pos, self.iln_strand = parse_description(
            sequence.description)

        logger.info(
            f"Processing SNP: {self.snp_id}, iln_snp: {self.iln_snp}, "
            f"iln_pos: {self.iln_pos}, iln_strand: {self.iln_strand} "
            f"probe_len: {self.probe_len}")

        self.got_manifest = True

    def filter_results(self, lenth_pct=60, ident_pct=97):
        # reset best hsp
        self.best_hit = None

        # hit represents a single database hit
        logger.debug(
            f"Got {len(self.queryresult.hits)} hits for {self.queryresult.id}")

        # filter results by score (query aligned)
        def filter_hsps(hsp: HSP):
            # score is a function of probe length
            if (hsp.score < self.probe_len * lenth_pct / 100 or
                    hsp.ident_pct < ident_pct):
                logger.debug(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Bad Score: {hsp.score} (ident_pct: {hsp.ident_pct})"
                )
                return False

            logger.debug(
                f"Keeping {hsp.hit_id}:{hsp.hit_range_all}: "
                f"Score: {hsp.score} (ident_pct: {hsp.ident_pct})"
            )

            return True

        filtered = self.queryresult.hsp_filter(filter_hsps)

        logger.info(f"Got {len(filtered.hits)} hits after filtering")

        if len(filtered.hits) == 0 or len(filtered.hsps) == 0:
            logger.warning(
                f"All alignments have been filtered out for "
                f"{self.queryresult.id}")

        elif len(filtered.hits) > 1 or len(filtered.hsps) > 1:
            for hsp in filtered.hsps:
                logger.debug(
                    f"{hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Score: {hsp.score} (ident_pct: {hsp.ident_pct})"
                )

            # try to sort and filter results
            filtered = self.sort_filtered(filtered)

            logger.warning(
                f"Got {len(filtered.hsps)} alignments after filtering for "
                f"{self.queryresult.id}")

        else:
            self.best_hit = filtered[0]
            self.status = "Found valid alignment after filtering"

        self.filtered = filtered
        self.is_filtered = True

    def sort_filtered(self, filtered):
        """Attempt to sort HSPs and get the best alignment"""

        # get first two HSPs
        hsp1, hsp2, *_ = sorted(
            filtered.hsps, key=lambda hsp: hsp.score, reverse=True)

        if hsp1.score >= hsp2.score:
            def filter_hsps(hsp: HSP):
                if hsp.score < hsp1.score:
                    logger.debug(
                        f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                        f"Bad Score: {hsp.score} (ident_pct: {hsp.ident_pct})"
                    )
                    return False

                logger.debug(
                    f"Keeping {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Score: {hsp.score} (ident_pct: {hsp.ident_pct})"
                )

                return True

            filtered = self.queryresult.hsp_filter(filter_hsps)

            logger.info(f"Got {len(filtered.hits)} hits after sorting")

        if len(filtered.hits) == 1 and len(filtered.hsps) == 1:
            self.best_hit = filtered[0]
            self.status = "Found valid alignment after filtering"

        else:
            self.status = "Too many alignments after filtering"

        return filtered

    def __discard_snp(self, message: str):
        """Returns a record for a discarded snp"""

        logger.error(message)
        logger.warning(f"Discarding {self}")
        self.status = message

        line = [
            self.snp_id, 0, 0, None, self.iln_snp, None,
            self.iln_strand, None, None, None]

        discarded = [
            self.snp_id,
            self.iln_snp,
            self.iln_strand,
            self.status]

        return line, discarded

    def _search_fragment(self, hsp: HSP):
        """Search for a SNP in HSP fragment"""

        found = None

        for fragment in hsp.fragments:
            start, end = fragment.query_range
            if self.iln_pos >= start and self.iln_pos <= end:
                logger.debug(f"found {self.iln_snp} in {fragment}")
                found = fragment

        if not found:
            raise BlatException(f"Cannot find '{self.iln_snp}' in fragment")

        return found

    def __process_hits(self, id2chromosome: dict):
        for i, hit in enumerate(self.filtered.hits):
            discarded = []
            alignments = []

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
                    f"Score {hsp.score}, ident_pct: {hsp.ident_pct}")

                # manage fragmented HSP
                if hsp.is_fragmented:
                    try:
                        hsp = self._search_fragment(hsp)
                    except BlatException as exc:
                        # process alignments
                        for fragment in hsp.fragments:
                            alignment = annotate_alignment(
                                fragment, chrom)
                            alignments.append(alignment)

                        line, discarded = self.__discard_snp(exc.value)
                        yield line, alignments, discarded
                        continue

                orient = check_strand(hsp)

                logger.info(f"Orient is '{orient}'")

                # the SNP position in the alignment, supposing no gap in
                # query sequence (mind to the query strand)
                if hsp.query_strand > 0:
                    snp_pos = self.iln_pos - hsp.query_start - 1
                else:
                    snp_pos = hsp.query_end - self.iln_pos

                # get and annotate alignment
                alignment = annotate_alignment(hsp, chrom)
                alignments.append(alignment)

                # check that SNP is in sequence
                if (snp_pos > hsp.query_end) or (snp_pos < hsp.query_start):
                    message = f"Can't find {self.iln_snp} in alignment"
                    line, discarded = self.__discard_snp(message)
                    yield line, alignments, discarded
                    continue

                # check that is letter is a N
                if alignment[0][snp_pos].upper() != 'N':
                    logger.error(alignment[:, snp_pos-5:snp_pos+6])
                    message = f"Cannot find SNP in position {snp_pos}"
                    line, discarded = self.__discard_snp(message)
                    yield line, alignments, discarded
                    continue

                if hsp.hit_strand > 0:
                    ref_pos = hsp.hit_start + snp_pos
                else:
                    ref_pos = hsp.hit_end - snp_pos - 1

                # this is 0-based index
                ref_allele = alignment[1][snp_pos].upper()

                # mind to reference strand
                if hsp.hit_strand < 0:
                    ref_allele = complement(ref_allele)

                logger.info(
                    f"Reference allele: {ref_allele} at "
                    f"{hit.id}:{ref_pos+1}"
                )

                try:
                    alt_allele = get_alt_allele(
                        self.iln_snp, ref_allele, hsp)

                except ValueError as exc:
                    logger.warning(
                        f"Cannot find alt_allele in illumina SNP: "
                        f"'{exc}")
                    message = "Allele doesn't match to reference"
                    line, discarded = self.__discard_snp(message)

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

                yield line, alignments, discarded

    def process_alignments(self, id2chromosome: dict):
        """Returns an output record and the processed alignment"""

        self.lines, self.alignments, discarded_snps = [], [], []

        if not self.is_filtered:
            raise BlatException(
                "You have to call 'filter_results()' before this method")

        if not self.got_manifest:
            raise BlatException(
                "You have to call 'read_sequence_manifest()' "
                "before this method")

        if not self.best_hit:
            message = "No valid alignments after filtering"
            line, discarded = self.__discard_snp(message)
            self.lines.append(line)
            discarded_snps.append(discarded)

        else:
            for line, alignments, discarded in self.__process_hits(
                    id2chromosome):
                self.lines.append(line)
                self.alignments += alignments

                if discarded:
                    discarded_snps.append(discarded)

        return self.lines, self.alignments, discarded_snps
