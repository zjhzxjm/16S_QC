"""
Author: xujm@realbio.cn
Ver:20160317

"""
import re


class SeqIndex():
    its_primer1_pattern = re.compile(r'TCCTCCGCTTATTGATATGC')
    its_primer2_pattern = re.compile(r'GCATCGATGAAGAACGCAGC')
    aflp_primer1_pattern = re.compile(r'GA[TC][GT][AG][GC][TG][CT][TA][AC][GC]AA[CT][GT][GC][TA]')
    aflp_primer2_pattern = re.compile(r'GA[TC][GT][AG][GC][TG][CT][TA][AC][GC]AA[CT][GT][GC][TA]')
    nifh_primer1_pattern = re.compile(r'TGCGA[CT]CC[GC]AA[AG]GC[GTC]GACTC')
    nifh_primer2_pattern = re.compile(r'AT[GC]GCCATCAT[CT]TC[AG]CCGGA')

    def __init__(self):
        pass

    out_barcode = {
        'realgene': 'ATCTCG',
    }

    primer = {
        'hiseq': {
            '16S': {
                'forward': 'GGACTACVVGGGTATCTAATC',
                'reverse': 'CCTACGGGRSGCAGCAG',
            },
            'ITS': {
                'forward': 'TCCTCCGCTTATTGATATGC',
                'reverse': 'GCATCGATGAAGAACGCAGC',
            },
            'ARCH': {
                'forward': 'GGACTACVVGGGTATCTAATC',
                'reverse': 'CCCTACGGGGYGCASCAG',
            },
            'AFLP': {
                'forward': 'GAYKRSKYWMSAAYKSW',
                'reverse': 'GAYKRSKYWMSAAYKSW',
            },
            'nifH': {
                'forward': 'TGCGAYCCSAARGCBGACTC',
                'reverse': 'ATSGCCATCATYTCRCCGGA',
            },
        },
        'miseq': {
            '16S': {
                'forward': 'GGACTACVVGGGTATCTAATC',
                'reverse': 'CCTACGGGRSGCAGCAG',
            },
            'ITS': {
                'forward': 'TCCTCCGCTTATTGATATGC',
                'reverse': 'GCATCGATGAAGAACGCAGC',
            },
            'ARCH': {
                'forward': 'GGACTACVVGGGTATCTAATC',
                'reverse': 'CCCTACGGGGYGCASCAG',
            },
        },
        'slim': {
            '16S': {
                'forward': 'TARGCCAAWACCKTACCA',
                'reverse': 'GAAACCTGGTTGATCCTG',
            },
        },
    }

    barcode = {
        'hiseq': ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT',
                  'GGCTAC', 'CTTGTA', 'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTAGAG', 'GTCCGC', 'GTGAAA', 'GTGGCC',
                  'GTTTCG', 'CGTACG', 'GAGTGG', 'GGTAGC', 'ACTGAT', 'ATGAGC', 'ATTCCT', 'CAAAAG', 'CAACTA', 'CACCGG'],
        'miseq': ['CCTAAA', 'TGCAGA', 'CCATCA', 'GTGGTA', 'ACTTTA', 'GAGCAA', 'TGTTGC', 'ATGTCC', 'AGGTAC', 'GTTACG',
                  'TACCGC', 'CGTAAG', 'ACAGCC', 'TGTCTC', 'GAGGAG', 'TACCGG', 'ATCTAG', 'CCAGGG', 'CACCTT', 'ATAGTT',
                  'GCACTT', 'TTAACT', 'CGCGGT', 'GAGACT']
    }
