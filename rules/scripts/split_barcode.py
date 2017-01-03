# -*- coding: utf-8 -*-
"""
Author: xujm@realbio.cn
Ver1.8
add nifH data type
Ver1.7
Support fuzzy search for barcode
Ver1.6
Add qibao AFLP data type
Ver1.5
Support repeat barcode between ITS and 16S
Ver1.4
Because hiseq ITS not synthesis again so it should recognized by primer
Ver1.3
Add out barcode and barcode name params
Ver1.2
1. According the seq description to judge out barcode exists or not
2. fix gzip error on python 3.5
Ver1.1
1. info.log only report barcode pair
Ver:1.0
init
"""

import argparse
import gzip
import logging
import os
import re
import subprocess

import fuzzysearch
from Bio import SeqIO

import setting

parser = argparse.ArgumentParser(description="Split samples from Illumina sequencing")
parser.add_argument('-a', '--fq1', type=str, dest='fq1', help='Read1 fastq file', required=True)
parser.add_argument('-b', '--fq2', type=str, dest='fq2', help='Read2 fastq file', required=True)
parser.add_argument('-n', '--barcodeName', type=str, dest='barcode_name', default="hiseq", help='Barcode name, \
                    the option for built-in: miseq and hiseq, default is hiseq. It will not affect when provide single \
                    barcode seq.')
parser.add_argument('-f', '--fuzzySearch', action='store_true', dest='fuzzy_search', help='Enable fuzzy search for \
                    inner barcode')
parser.add_argument('-s', '--sampleConfig', type=str, dest='sample_config', help='Sample barcode configuration info',
                    required=True)
parser.add_argument('-o', '--outBarcode', type=str, dest='out_barcode', help='Default is built-in \
                    realgene seq')
parser.add_argument('-w', '--workDir', type=str, dest='work_dir', default=".", help='Work directory, default is ./')
parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help='Enable debug info')
parser.add_argument('--version', action='version', version='1.8')


class RawFastqPairInfo:
    def __init__(self, ob_read1, ob_read2, outbarcode, barcodename, fuzzsearch=False):
        self.ob_read1 = ob_read1
        self.ob_read2 = ob_read2
        self.out_barcode = outbarcode
        self.barcode_name = barcodename
        self.fuzzy_search = fuzzsearch

    def get_barcode_pair(self):
        fq1_seq = str(self.ob_read1.seq)
        fq2_seq = str(self.ob_read2.seq)

        data_type = ''
        if self.barcode_name == "hiseq":
            fq1_bar = fq1_seq[:6]
            fq2_bar = fq2_seq[:6]
            barcode = self.barcode_name
            if re.search(setting.SeqIndex.its_primer1_pattern, fq1_seq) and re.search(
                    setting.SeqIndex.its_primer2_pattern, fq2_seq):
                fq1_bar = fq2_seq[:6]
                fq2_bar = fq1_seq[:6]
                barcode = 'miseq'
                data_type = 'ITS-'
            elif re.search(setting.SeqIndex.aflp_primer1_pattern, fq1_seq) and re.search(
                    setting.SeqIndex.aflp_primer2_pattern, fq2_seq):
                fq1_bar = fq1_seq[:6]
                fq2_bar = fq2_seq[:6]
                barcode = 'hiseq'
                data_type = 'AFLP-'
            elif re.search(setting.SeqIndex.nifh_primer1_pattern, fq1_seq) and re.search(
                    setting.SeqIndex.nifh_primer2_pattern, fq2_seq):
                fq1_bar = fq1_seq[:6]
                fq2_bar = fq2_seq[:6]
                barcode = 'hiseq'
                data_type = 'nifH-'

        elif self.barcode_name == "miseq":
            fq1_bar = fq2_seq[:6]
            fq2_bar = fq1_seq[:6]
            barcode = self.barcode_name
        else:
            logging.critical('barcode {0} not exists'.format(self.barcode_name))

        try:
            f_barcode = "F" + str(setting.SeqIndex.barcode[barcode].index(fq1_bar) + 1)
        except ValueError:
            pos = self.fuzzy_search_barcode(barcode, fq1_bar)
            if self.fuzzy_search and pos:
                f_barcode = "F" + pos
            else:
                f_barcode = ''

        try:
            r_barcode = "R" + str(setting.SeqIndex.barcode[barcode].index(fq2_bar) + 1)
        except ValueError:
            pos = self.fuzzy_search_barcode(barcode, fq2_bar)
            if self.fuzzy_search and pos:
                r_barcode = "R" + pos
            else:
                r_barcode = ''

        if f_barcode and r_barcode:
            return data_type + f_barcode + "+" + r_barcode

    def fuzzy_search_barcode(self, barcodes_name, inner_barcode):
        """
        Fuzzy search inner barcode in barcodes table, for sequencing company's problem
        :param barcodes_name:
        :param inner_barcode:
        :return: the position of barcodes table
        """
        position = 0
        flag = 0
        tag = ''
        while position < len(setting.SeqIndex.barcode[barcodes_name]):
            if flag > 1:
                tag = ''
                break
            if fuzzysearch.find_near_matches(setting.SeqIndex.barcode[barcodes_name][position], inner_barcode, 1, 1, 1, 1):
                tag = str(position + 1)
                flag += 1
            position += 1
        return tag

    def is_need_out_barcode(self):
        """
        Judge out barcode in each reads pair
        Args:
        Returns: True or False

        """
        try:
            if re.match(r'[ATCG]', self.ob_read1.description.split()[1].split(":")[-1].strip()):
                if fuzzysearch.find_near_matches(self.out_barcode, self.ob_read1.description, 1, 0, 0, 1) \
                        and fuzzysearch.find_near_matches(self.out_barcode, self.ob_read2.description, 1, 0, 0, 1):
                    return True
                else:
                    return False
            else:
                return True
        except IndexError:
            return True


class RawFastqSingleInfo:
    def __init__(self, ob_read):
        self.ob_read = ob_read

    def get_barcode(self):
        fq_seq = str(self.ob_read.seq)
        fq_bar = fq_seq[:12]
        return fq_bar


class Sample:
    def __init__(self, sam_barcode, work_dir='.'):
        d_dir = {}
        l_check_barcode_type = []
        with open(sam_barcode) as f:
            for line in f:
                (project, sample, barcode, data_type, lib_type) = line.strip().split()[:5]
                code = subprocess.call(['mkdir', '-p', work_dir + "/" + project + "/" + sample + "_" + data_type])
                if code:
                    logging.error("Can't make filefoder: %s/%s/%s" % (work_dir, project, sample))
                if data_type == 'ITS':
                    d_dir['ITS-' + barcode] = work_dir + "/" + project + "/" + sample + "_" + data_type
                elif data_type == 'AFLP':
                    d_dir['AFLP-' + barcode] = work_dir + "/" + project + "/" + sample + "_" + data_type
                elif data_type == 'nifH':
                    d_dir['nifH-' + barcode] = work_dir + "/" + project + "/" + sample + "_" + data_type
                else:
                    d_dir[barcode] = work_dir + "/" + project + "/" + sample + "_" + data_type
                l_check_barcode_type.append(len(barcode.split("+")))

        if len(list(set(l_check_barcode_type))) == 1 and list(set(l_check_barcode_type))[0] == 2:
            barcode_type = "pair"
        elif len(list(set(l_check_barcode_type))) == 1 and list(set(l_check_barcode_type))[0] == 1:
            barcode_type = "single"
        else:
            print("ERROR:illegal barcode in sam_barcode.all")

        subprocess.call(['mkdir', '-p', work_dir + '/Unalign'])
        self.d_dir = d_dir
        self.barcode_type = barcode_type


if __name__ == '__main__':
    args = parser.parse_args()
    fq1 = os.path.abspath(args.fq1)
    fq2 = os.path.abspath(args.fq2)
    fq1_name = fq1.split('/')[-1]
    fq2_name = fq2.split('/')[-1]
    barcode_name = args.barcode_name
    if args.fuzzy_search:
        fuzzy_search = True
    else:
        fuzzy_search = False
    work_path = os.path.abspath(args.work_dir)

    if args.out_barcode:
        out_barcode = args.out_barcode
    else:
        out_barcode = setting.SeqIndex.out_barcode['realgene']

    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG,
            format="[%(asctime)s]%(name)s:%(levelname)s:%(message)s",
            filename=work_path + "/debug.log"
        )
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s]%(name)s:%(levelname)s:%(message)s",
            filename=work_path + "/" + fq1_name.replace("R1.", "") + ".info.log"
        )
    logging.debug("Start running")

    class_sample = Sample(args.sample_config, work_path)

    if re.findall(r'gz', fq1):
        F_fq1 = gzip.open(fq1, "rt")
    else:
        F_fq1 = open(fq1)
    if re.findall(r'gz', fq2):
        F_fq2 = gzip.open(fq2, "rt")
    else:
        F_fq2 = open(fq2)

    O_fq1 = {}
    O_fq2 = {}
    for (k, v) in class_sample.d_dir.items():
        O_fq1[k] = open(v + "/" + fq1_name + ".filtered", "w")
        O_fq2[k] = open(v + "/" + fq2_name + ".filtered", "w")
    O_unalign_fq1 = open(work_path + "/Unalign/" + fq1_name + ".unalign", "w")
    O_unalign_fq2 = open(work_path + "/Unalign/" + fq2_name + ".unalign", "w")

    fq1_iter = SeqIO.parse(F_fq1, "fastq")
    fq2_iter = SeqIO.parse(F_fq2, "fastq")

    d_count = {'total': 0, 'out_total': 0}
    while True:
        try:
            record_fq1 = next(fq1_iter)
            record_fq2 = next(fq2_iter)
            if class_sample.barcode_type == "pair":
                # logging.debug("Pair run lib_type {0}".format(class_sample.lib_type))
                class_fastq_pair = RawFastqPairInfo(record_fq1, record_fq2, out_barcode, barcode_name,
                                                    fuzzsearch=fuzzy_search)

                d_count['total'] += 1
                if class_fastq_pair.is_need_out_barcode():
                    # Fetch our out barcode
                    d_count['out_total'] += 1
                    barcode_pair = class_fastq_pair.get_barcode_pair()
                    if barcode_pair:
                        try:
                            # logging.debug("Our seq %s" % class_sample.d_dir[barcode_pair])
                            O_fq1[barcode_pair].write(record_fq1.format("fastq"))
                            O_fq2[barcode_pair].write(record_fq2.format("fastq"))
                            # logging.debug(record_fq1.format("fastq"))
                            # try:
                            #     d_count[class_sample.d_dir[barcode_pair]] += 1
                            # except KeyError:
                            #     d_count[class_sample.d_dir[barcode_pair]] = 1
                        except KeyError:
                            O_unalign_fq1.write(record_fq1.format("fastq"))
                            O_unalign_fq2.write(record_fq2.format("fastq"))
                            try:
                                d_count[barcode_pair] += 1
                            except KeyError:
                                d_count[barcode_pair] = 1
                    else:
                        O_unalign_fq1.write(record_fq1.format("fastq"))
                        O_unalign_fq2.write(record_fq2.format("fastq"))
            elif class_sample.barcode_type == "single":
                logging.debug("Pair run lib_type {0}".format(class_sample.lib_type))
                class_fastq_single = RawFastqSingleInfo(record_fq2)
                barcode_single = class_fastq_single.get_barcode()
                try:
                    O_fq1[barcode_single].write(record_fq1.format("fastq"))
                    O_fq2[barcode_single].write(record_fq2.format("fastq"))
                except KeyError:
                    O_unalign_fq1.write(record_fq1.format("fastq"))
                    O_unalign_fq2.write(record_fq2.format("fastq"))
            else:
                logging.error("Null barcode type from sub function")
                break
        except StopIteration:
            break

    for k in class_sample.d_dir.keys():
        O_fq1[k].close()
        O_fq2[k].close()

    if class_sample.barcode_type == "pair":
        for (k, v) in d_count.items():
            if k == "total":
                logging.debug("The total reads number: %d" % v)
            elif k == "out_total":
                logging.debug("Total reads number with our out barcode: %d" % v)
            else:
                logging.info("%s: %d" % (k, v))

    logging.debug('End running')
