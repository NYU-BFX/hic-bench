#! /usr/bin/env python

"""
extract_reads.py
Created by Tim Stuart
"""

import pysam


def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n


def extract_reads(options):
    n = get_names(options.names)
    bamfile = pysam.AlignmentFile(options.bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()
    out = pysam.Samfile(options.out, 'wb', header=header)
    for name in n:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Extract reads by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='list of read names to extract', required=True)
    parser.add_argument('-o', '--out', help='file name for extracted alignments', required=True)   
    options = parser.parse_args()
    extract_reads(options)
