"""
Author: xujm@realbio.cn
Ver:20160317

"""


class SeqIndex():
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
        },
        'miseq': {
            '16S': {
                'forward': '',
                'reverse': '',
            },
            'ITS': {
                'forward': '',
                'reverse': '',
            },
            'ARCH': {
                'forward': '',
                'reverse': '',
            },
        }
    }

    barcode = {
        'realgene': ['ATCACG', 'CGATGT', 'TTAGGC', 'TGACCA', 'ACAGTG', 'GCCAAT', 'CAGATC', 'ACTTGA', 'GATCAG', 'TAGCTT',
                     'GGCTAC', 'CTTGTA', 'AGTCAA', 'AGTTCC', 'ATGTCA', 'CCGTCC', 'GTAGAG', 'GTCCGC', 'GTGAAA', 'GTGGCC',
                     'GTTTCG', 'CGTACG', 'GAGTGG', 'GGTAGC', 'ACTGAT', 'ATGAGC', 'ATTCCT', 'CAAAAG', 'CAACTA', 'CACCGG']
    }
