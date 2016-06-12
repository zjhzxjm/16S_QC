"""
Author: xujm@realbio.cn
Ver:

"""
# -*- coding: utf-8 -*- \#

import os, re, sys
import argparse
import logging

parser = argparse.ArgumentParser(description="")
parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help='Enable debug info')

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
