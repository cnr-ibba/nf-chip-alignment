#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 16:39:44 2023

@author: Paolo Cozzi <bunop@libero.it>
"""

import re
import logging

from Bio.SearchIO._model.query import QueryResult
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import SNP2BASES
from .utils import complement

# Get an instance of a logger
logger = logging.getLogger(__name__)

manifest_pattern = re.compile(
    "iln_snp: (.*), iln_pos: (.*), iln_strand: (.*)")


def parse_description(description):
    match = re.search(manifest_pattern, description)
    iln_snp, iln_pos, iln_strand = match.groups()
    return iln_snp, int(iln_pos), iln_strand


def rename_hit(hit):
    # try to determine chromosome name
    match_chrom = re.search(r"chromosome (\w*),", hit.description)
    match_scaff = re.search(r"(scaffold_?[0-9]+),", hit.description)
    match_unknw = re.search(r"(unplaced_?[0-9]+),", hit.description)

    if match_chrom:
        hit.id = match_chrom.groups()[0]

    elif match_scaff:
        hit.id = match_scaff.groups()[0]

    elif match_unknw:
        hit.id = match_unknw.groups()[0]

    return hit


def check_strand(hsp):
    """Return orient respect to the hit_strand"""

    orient = hsp.query_strand * hsp.hit_strand

    if orient > 0:
        return "forward"

    else:
        return "reverse"


def get_alt_allele(snp, ref, hsp):
    """Get ALT allele giving the ref and orientation"""

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


class BlastResult():
    queryresult = None
    snp_id = None
    iln_snp = None
    iln_pos = None
    iln_strand = None
    probe_len = None
    status = None

    def __init__(self, queryresult: QueryResult = None):
        if queryresult is not None:
            self._parse_queryresult(queryresult)

    def _parse_queryresult(self, queryresult):
        self.snp_id = queryresult.id
        self.iln_snp, self.iln_pos, self.iln_strand = parse_description(
            queryresult.description
        )
        self.probe_len = queryresult.seq_len

        # parse chromosome names and change ids
        if queryresult.hits:
            queryresult = queryresult.hit_map(rename_hit)

        self.queryresult = queryresult

    def __repr__(self):
        return (
            f"BlastResult: {self.snp_id}, iln_snp: {self.iln_snp}, "
            f"iln_pos: {self.iln_pos}, iln_strand: {self.iln_strand}, "
            f"probe_len: {self.probe_len}")

    def filter_results(self, length_pct=95, ident_pct=97):
        # reset best hsp
        self.best_hit = None

        # hit represents a single database hit
        logger.info(
            f"Got {len(self.queryresult.hits)} hits for {self.queryresult.id}")

        # filter results by alignment length and identity
        def filter_hsps(hsp):
            if hsp.is_fragmented:
                logger.warning(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Found {len(hsp.fragments)} fragments"
                )
                return False

            query_length_pct = hsp.query_span / self.probe_len * 100
            query_ident_pct = hsp.ident_num / hsp.query_span * 100

            if query_length_pct < length_pct or query_ident_pct < ident_pct:
                logger.debug(
                    f"Filtering out {hsp.hit_id}:{hsp.hit_range_all}: "
                    f"Bad Score: {hsp.bitscore} (ident_pct: {query_ident_pct},"
                    f" length_pct: {query_length_pct})"
                )
                return False

            logger.debug(
                f"Keeping {hsp.hit_id}:{hsp.hit_range_all}: "
                f"Score: {hsp.bitscore} (ident_pct: {query_ident_pct},"
                f" length_pct: {query_length_pct})"
            )

            return True

        # apply first filter on length
        filtered = self.queryresult.hsp_filter(filter_hsps)

        # filter by HSP score:
        # get best score and filter all HSPs with the same score
        if len(filtered) > 0:
            evalues = [hsp.evalue for hsp in filtered.hsps]
            best_score = min(evalues)

        # apply second filter on HSP score
        filtered = filtered.hsp_filter(lambda hsp: hsp.evalue <= best_score)

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
                    f"Score: {hsp.bitscore} (ident_pct: "
                    f"{hsp.ident_num / self.probe_len * 100})"
                )

            self.status = "Too many alignments after filtering"

            logger.warning(
                f"Got {len(filtered.hsps)} alignments after filtering for "
                f"{self.queryresult.id}")

        else:
            self.best_hit = filtered[0]
            self.status = "Found valid alignment after filtering"

        self.filtered = filtered
        self.is_filtered = True

    def __discard_snp(self, status=None):
        """Returns a record for a discarded snp"""

        if status:
            logger.error(status)
            logger.warning(f"Discarding {self}")
            self.status = status

        line = [
            self.snp_id, 0, 0, None, self.iln_snp, None,
            self.iln_strand, None, None, None]

        discarded = [
            self.snp_id,
            self.iln_snp,
            self.iln_strand,
            self.status]

        return line, discarded

    def __find_snp_pos(self, hsp):
        """Search SNP in alignment query. Return 0-based position"""

        # the SNP position in the alignment (mind to the query strand)
        if hsp.query_strand > 0:
            snp_pos = self.iln_pos - hsp.query_start - 1

            # test for gapped probes
            idx = hsp.query.seq.find(SNP2BASES[self.iln_snp])
            gaps = hsp.query.seq.count("-", 0, idx)

            # adjust snp_pos
            snp_pos += gaps

        else:
            snp_pos = hsp.query_end - self.iln_pos

            # test for gapped probes
            idx = hsp.query.seq.find(SNP2BASES[self.iln_snp])
            gaps = hsp.query.seq.count("-", idx)

            # adjust snp_pos
            snp_pos -= gaps

        return snp_pos

    def __find_ref_pos(self, hsp, snp_pos):
        """Search for REF position on HIT sequence"""

        if hsp.hit_strand > 0:
            ref_pos = hsp.hit_start + snp_pos

            # test for gapped reference
            gaps = hsp.hit.seq.count("-", 0, snp_pos)

            # adjust snp_pos
            ref_pos -= gaps

        else:
            ref_pos = hsp.hit_end - snp_pos - 1

            # test for gapped reference
            gaps = hsp.hit.seq.count("-", snp_pos)

            # adjust snp_pos
            ref_pos += gaps

        return ref_pos

    def __reject_adjacent_gaps(self, hsp):
        """check and discard SNP with adjacent gaps when reference sequence
        is compatible with the SNP call"""

        # the SNP position in the alignment (mind to the query strand)
        idx = hsp.query.seq.find(SNP2BASES[self.iln_snp])

        # check and discard SNP after/before gaps, with reference
        # sequences compatible with self.iln_snp: in those cases
        # I can't determine a unique position. See this example
        # GGAA-YCAGTT 250506CS3900371000001_1255.1
        # GGAACTCAGTT 11:35291071-35291193
        # ||||  |||||

        # test for gap before SNP
        if hsp.query.seq[idx-1] == '-':
            base = hsp.hit.seq[idx-1]
            if base in self.iln_snp.split("/"):
                return True

        # test for gap after SNP
        if hsp.query.seq[idx+1] == '-':
            base = hsp.hit.seq[idx+1]
            if base in self.iln_snp.split("/"):
                return True

        # check and discard SNP after/before gaps, with reference
        # allele compatible with SNP neighbour. Something like
        # TAAGYCGAGG
        # TAAGC-GAGG
        # ||||  ||||

        # ensure no cases with gaps on reference sequence
        # in this case
        if hsp.hit.seq[idx-1] == '-':
            base = hsp.hit.seq[idx]
            if base == hsp.query.seq[idx-1]:
                return True

        if hsp.hit.seq[idx+1] == '-':
            base = hsp.hit.seq[idx]
            if base == hsp.query.seq[idx+1]:
                return True

        return False

    def __process_hits(self):
        discarded = []

        for i, hit in enumerate(self.filtered.hits):
            logger.info(f"Processing hit {i}: {hit.id} for {self.snp_id}")

            # hsp represents region(s) of significant alignments between
            # query and hit sequences
            logger.debug(f"Got {len(hit.hsps)} hsp for {hit.id}")

            for j, hsp in enumerate(hit.hsps):
                logger.info(
                    f"Hsp {j}: has {len(hsp.fragments)} fragments. "
                    f"Query range {hsp.query_range_all} "
                    f"({hsp.query_strand_all}), "
                    f"Hit range {hsp.hit_range_all} ({hsp.hit_strand_all}), "
                    f"Score {hsp.bitscore}, ident_pct "
                    f"{hsp.ident_num / self.probe_len * 100}")

                orient = check_strand(hsp)

                logger.info(f"Orient is '{orient}'")

                snp_pos = self.__find_snp_pos(hsp)

                # get and annotate alignment
                alignment = hsp.aln
                alignment[-1].id = "{0}:{1}-{2}".format(
                    hit.id,
                    hsp.hit_range[0],
                    hsp.hit_range[1])

                # add the similarity annotation
                annotation = SeqRecord(
                    Seq(hsp.aln_annotation['similarity']), id='')
                alignment.append(annotation)

                # check that SNP is in sequence (mind that snp_pos is scaled
                # by query start)
                if ((snp_pos >= (hsp.query_end - hsp.query_start)) or
                        (snp_pos < hsp.query_start)):
                    message = f"Can't find SNP in alignment"
                    line, discarded = self.__discard_snp(message)
                    yield line, alignment, discarded
                    continue

                # check snp position with IUPAC ambiguity codes
                if alignment[0][snp_pos].upper() != SNP2BASES[self.iln_snp]:
                    message = f"Cannot find the expected SNP in sequence"
                    logger.error(alignment[:, snp_pos-5:snp_pos+6])
                    line, discarded = self.__discard_snp(message)
                    yield line, alignment, discarded
                    continue

                # check and discard SNP after/before gaps, with reference
                # sequences compatible with self.iln_snp: in those cases
                # I can't determine a unique position.
                if self.__reject_adjacent_gaps(hsp):
                    message = (
                        f"Cannot determine a unique SNP position")
                    logger.error(alignment[:, snp_pos-5:snp_pos+6])
                    line, discarded = self.__discard_snp(message)
                    yield line, alignment, discarded
                    continue

                ref_pos = self.__find_ref_pos(hsp, snp_pos)

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

                except ValueError:
                    logger.warning(
                        f"Cannot find alt_allele in illumina SNP: "
                        f"'{ref_allele}' not in {self.iln_snp}")
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
                        self.snp_id, hit.id, ref_pos+1, alleles, self.iln_snp,
                        get_illumina_forward(self.iln_snp, orient),
                        self.iln_strand, orient, ref_allele, alt_allele]

                yield line, alignment, discarded

    def process_alignments(self):
        """Returns an output record and the processed alignment"""

        self.lines, self.alignments, discarded_snps = [], [], []

        if self.is_filtered and not self.best_hit:
            logger.warning(f"Discarding {self}")
            # return an empty alignment
            line, discarded = self.__discard_snp()
            self.lines.append(line)
            discarded_snps.append(discarded)

        else:
            for line, alignment, discarded in self.__process_hits():
                self.lines.append(line)
                self.alignments.append(alignment)

                if discarded:
                    discarded_snps.append(discarded)

        return self.lines, self.alignments, discarded_snps
