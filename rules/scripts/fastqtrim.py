"""
Author: xujm@realbio.cn
Ver:1.0
init
"""
# -*- coding: utf-8 -*- \#

import os, re, sys
import argparse
import logging
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument('-i', '--input', type=str, dest='input', help='fastq file', required=True)
parser.add_argument('-o', '--output', type=str, dest='output', help='filtered fastq out', required=True)
parser.add_argument('-b', '--bases', type=int, dest='bases', default=50, help='the number of bases, default is 50')
parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help='Enable debug info')
parser.add_argument('--version', action='version', version='1.0')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="[%(asctime)s]%(name)s:%(levelname)s:%(message)s",
            filename='debug.log'
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s]%(name)s:%(levelname)s:%(message)s",
            filename='info.log'
        )

if __name__ == '__main__':
    args = parser.parse_args()
    fq = args.input
    out = args.output
    bases = args.bases

    fq_iter = SeqIO.parse(open(fq), "fastq")
    O_fq = open(out, "w")
    trimmed_reads = (rec[:-bases] for rec in fq_iter)
    SeqIO.write(trimmed_reads, out, "fastq")
