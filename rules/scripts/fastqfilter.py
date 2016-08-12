"""
Author: xujm@realbio.cn
Ver:1.2
Fastqfilter sample name error bug #9
Ver:1.1
1。fix output relative path error
2. remove info.log output
Ver:1.0
init
"""
# -*- coding: utf-8 -*- \#

import os, re, sys
import argparse
import logging
from Bio import SeqIO
from numpy import mean

parser = argparse.ArgumentParser(description="")
parser.add_argument('-i', '--input', type=str, dest='input', help='fastq file', required=True)
parser.add_argument('-o', '--output', type=str, dest='output', help='filtered fastq out', required=True)
parser.add_argument('-q', '--qmin', type=int, dest='qmin', default=20, help='min quality score, default is 20')
parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help='Enable debug info')
parser.add_argument('--version', action='version', version='1.1')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="[%(asctime)s]%(name)s:%(levelname)s:%(message)s",
            filename='debug.log'
        )

    fq = args.input
    out = os.path.abspath(args.output)
    qmin = args.qmin

    recover_sample_name = lambda x: x.replace('_AFLP', '').replace('_ITS', '').replace('_16S', '')
    sample_name = " " + recover_sample_name(out.split("/")[-2])
    fq_iter = SeqIO.parse(open(fq), "fastq")
    O_fq = open(out, "w")
    count = 0
    while True:
        try:
            record = next(fq_iter)
            logging.debug(type(mean(record.letter_annotations["phred_quality"])))
            if qmin > mean(record.letter_annotations["phred_quality"]):
                continue
            count += 1
            record.description += sample_name + "_" + str(count)
            O_fq.write(record.format("fastq"))
        except StopIteration:
            break
    O_fq.close()
