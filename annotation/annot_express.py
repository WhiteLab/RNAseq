#!/usr/bin/env python3

import sys


def make_index(index_ref):
    ind = {}
    index = open(index_ref, 'r')
    next(index)
    for line in index:
        ann = line.rstrip('\n').split('\t')
        ind[ann[0]] = {}
        ind[ann[0]]['name'] = ann[1]
        ind[ann[0]]['type'] = ann[2]
    index.close()
    return ind


def annot_express(index_ref, sample):
    ind = make_index(index_ref)
    table = open(sample + '.express_quantification.txt', 'r')
    head = next(table)
    out_fn = sample + '.express_quantification_annotated.txt'
    out_fh = open(out_fn, 'w')
    out_fh.write('name\ttype\t' + head)
    for line in table:
        info = line.split('\t')
        if float(info[3]) > 0:
            out_fh.write(ind[info[1]]['name'] + '\t ' + ind[info[1]]['type'] + '\t' + line)
        else:
            sys.stderr.write('Skipped ' + ind[info[1]]['name'] + ' ' + ind[info[1]]['type'] + ' ' + info[1]
                             + ' no reads!\n')
    table.close()
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-i', '--index', action='store', dest='index', help='Reference file with RNA gene names,'
                                                                                ' type and trasncript ids')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample name prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    index = inputs.index
    sample = inputs.sample

    annot_express(index, sample)

