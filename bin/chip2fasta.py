#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 12:33:06 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

import logging
import argparse

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
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    for i, record in enumerate(read_Manifest(args.input, delimiter=",")):
        print(record)
