#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 13:01:27 2022

@author: Paolo Cozzi <bunop@libero.it>
"""

# https://droog.gs.washington.edu/parc/images/iupac.html
IUPAC_AMBIGUITY_CODES = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'T',
    'M': 'A/C',
    'R': 'A/G',
    'W': 'A/T',
    'S': 'C/G',
    'Y': 'C/T',
    'K': 'G/T',
    'V': 'A/C/G',
    'H': 'A/C/T',
    'D': 'A/G/T',
    'B': 'C/G/T',
    'N': 'G/A/T/C'
}

SNP2BASES = {
    'A/C': 'M',
    'C/A': 'M',
    'A/G': 'R',
    'G/A': 'R',
    'A/T': 'W',
    'T/A': 'W',
    'C/G': 'S',
    'G/C': 'S',
    'C/T': 'Y',
    'T/C': 'Y',
    'G/T': 'K',
    'T/G': 'K'
}
