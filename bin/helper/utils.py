#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:13:51 2021

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
"""

import io
import re
import gzip
import logging
import pathlib

# Get an instance of a logger
logger = logging.getLogger(__name__)


def sanitize(
        word: str,
        chars=['.', ",", "-", "/", "#"],
        check_mongoengine=True) -> str:
    """Sanitize a word by removing unwanted characters and lowercase it.

    Args:
        word (str): the word to sanitize
        chars (list): a list of characters to remove
        check_mongoengine (bool): true to add '_' after a mongoengine reserved
            word

    Returns:
        str: the sanitized word
    """

    # remove unwanted characters from word by putting spaces
    pattern = "".join(chars)
    tmp = re.sub(r'[%s]' % (pattern), ' ', word)

    # remove spaces from column name and lowercase all
    sanitized = re.sub(r"\s+", "_", tmp).lower()

    if check_mongoengine:
        if sanitized in ['size', 'type']:
            sanitized += "_"

    # remove starting sanitized char (can't be used with namedtuple)
    if sanitized.startswith("_"):
        sanitized = sanitized[1:]

    return sanitized


def text_or_gzip_open(path: str, mode=None):
    if pathlib.Path(path).suffix == '.gz':
        if not mode:
            mode = 'rt'

        logger.debug(f"Gzip detected for {path}")
        return gzip.open(path, mode=mode)

    else:
        if not mode:
            mode = 'r'

        return open(path, mode=mode)


def complement(genotype: str):
    bases = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "/": "/",
        "-": "-"
    }

    result = ""

    for base in genotype:
        result += bases[base]

    return result


def skip_comments(handle: io.TextIOWrapper, comment_char="#") -> (int, list):
    """
    Ignore comments lines from a open file handle. Return the stream position
    immediately after the comments and all the comment lines in a list.

    Parameters
    ----------
    handle : io.TextIOWrapper
        An open file handle.
    comment_char : TYPE, optional
        The comment character used in file. The default is "#".

    Returns
    -------
    (int, list)
        The stream position after the comments and the ignored lines as a list.
    """

    # track skipped lines
    skipped = list()

    # read first line
    line = handle.readline().strip()

    # search for comments in file
    while line[0] == "#":
        logger.debug(f"Skipping: {line}")
        skipped.append(line)
        position = handle.tell()

        # read another line
        line = handle.readline().strip()

    # the position returned is the one before the one I want
    return position, skipped
