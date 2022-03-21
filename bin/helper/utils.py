#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 14:13:51 2021

@author: Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
"""

import gzip
import logging
import pathlib

# Get an instance of a logger
logger = logging.getLogger(__name__)


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
        "/": "/"
    }

    result = ""

    for base in genotype:
        result += bases[base]

    return result
