"""
Author: xujm@realbio.cn
Ver:

"""
# -*- coding: utf-8 -*- \#

import os, re, sys
import argparse
import logging
from Bio import SeqIO

parser = argparse.ArgumentParser(description="")
parser.add_argument('-a', '--fq1', type=str, dest='fq1', help='Read1 fastq file', required=True)
parser.add_argument('-b', '--fq2', type=str, dest='fq2', help='Read2 fastq file', required=True)
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
            filename="info.log"
        )

if __name__ == '__main__':
    args = parser.parse_args()
    fq1 = os.path.abspath(args.fq1)
    fq2 = os.path.abspath(args.fq2)

    out1 = fq1 + ".format"
    out2 = fq2 + ".format"

    N_max = 3
    len_min = 250

    O_fq1 = open(out1, "w")
    O_fq2 = open(out2, "w")

    d_name = {}

    for rec in SeqIO.parse(open(fq1), "fastq"):
        if len(rec.seq) >= len_min and len(re.findall("N", str(rec.seq))) < N_max:
            d_name[rec.name] = 0

    for rec in SeqIO.parse(open(fq2), "fastq"):
        try:
            logging.debug("fq2 common name:" + rec.name)
            logging.debug("fq2 len:" + str(len(rec.seq)))
            logging.debug("fq2 N num:" + str(len(re.findall("N", str(rec.seq)))))
            if len(rec.seq) < len_min or len(re.findall("N", str(rec.seq))) >= N_max:
                continue
            d_name[rec.name] += 1
            O_fq2.write(rec.format("fastq"))
        except KeyError:
            continue

    for rec in SeqIO.parse(open(fq1), "fastq"):
        logging.debug("fq1 len:" + str(len(rec.seq)))
        logging.debug("fq1 N num:" + str(len(re.findall("N", str(rec.seq)))))
        try:
            if d_name[rec.name]:
                logging.debug("fq1 name flag:" + str(d_name[rec.name]))
                O_fq1.write(rec.format("fastq"))
        except KeyError:
            continue

    O_fq1.close()
    O_fq2.close()
