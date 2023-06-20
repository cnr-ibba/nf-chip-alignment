#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 16:39:57 2021

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
"""

import csv
import logging
import collections
import datetime

from pathlib import Path
from typing import Union
from dateutil.parser import parse as parse_date

from .utils import (
    sanitize, text_or_gzip_open, skip_comments)

# Get an instance of a logger
logger = logging.getLogger(__name__)


def _search_in_header(header: list, term: str) -> list:
    return list(
        filter(
            lambda record: record.startswith(term),
            header
        )
    )


def search_manifactured_date(header: list) -> Union[datetime.datetime, None]:
    """Grep manifactured date from affymetrix header

    Args:
        header (list): affymetrix header section

    Returns:
        datetime.datetime: a datetime object
    """

    records = _search_in_header(header, "#%create_date=")
    date = None

    if records:
        record = records[0].split("=")
        date = parse_date(record[-1], fuzzy=True)

    return date


def search_n_samples(header: list) -> int:
    """Grep number of samples in affymetrix reportfile

    Args:
        header (list): affymetrix header section

    Returns:
        int: the number of samples in file
    """

    records = _search_in_header(header, "##samples-per-snp=")
    n_samples = None

    if records:
        n_samples = records[0].split("=")[1]
        n_samples = int(n_samples)

    return n_samples


def search_n_snps(header: list) -> int:
    """Grep number of SNPs in affymetrix reportfile

    Args:
        header (list): affymetrix header section

    Returns:
        int: the number of SNPs in file
    """

    records = _search_in_header(header, "##snp-count=")
    n_snp = None

    if records:
        n_snp = records[0].split("=")[1]
        n_snp = int(n_snp)

    return n_snp


def read_Manifest(path: Path, delimiter: str = ",") -> collections.namedtuple:
    """
    Open an affymetrix manifest file and yields records as namedtuple. Add
    an additional column for manifacured date (when SNP is recorded in
    datafile)

    Parameters
    ----------
    path : Path
        The position of manifest file.
    delimiter : str, optional
        field delimiter. The default is ",".

    Yields
    ------
    record : collections.namedtuple
        A single SNP record from manifest.
    """

    with text_or_gzip_open(path) as handle:
        position, skipped = skip_comments(handle)

        # go back to header section
        handle.seek(position)

        # now read csv file
        reader = csv.reader(handle, delimiter=delimiter)

        # get header
        header = next(reader)

        # sanitize column names
        header = [sanitize(column) for column in header]

        # replace sequence column name like illumina does
        columns = {
            'flank': 'sourceseq',
            'affy_snp_id': 'name'
        }

        for key, value in columns.items():
            if key in header:
                idx = header.index(key)
                header[idx] = value

        # ok try to get the manifatcured date
        date = search_manifactured_date(skipped)

        # add date to header
        header.append("date")

        logger.info(header)

        # define a datatype for my data
        SnpChip = collections.namedtuple("SnpChip", header)

        # add records to data
        for record in reader:
            # add date to record
            record.append(date)

            # set null values to items
            record = [col if col != '---' else None for col in record]

            # mind to null position
            if not record[header.index("chromosome")]:
                record[header.index("chromosome")] = "0"

            if not record[header.index("physical_position")]:
                record[header.index("physical_position")] = 0

            # convert into collection
            record = SnpChip._make(record)
            yield record
