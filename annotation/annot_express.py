#!/usr/bin/env python

import sys

index = open(sys.argv[1], 'r')

table = open(sys.argv[2], 'r')

ind = {}

head = next(index)
for line in index:
    ann = line.rstrip('\n').split('\t')
    ind[ann[0]] = {}
    ind[ann[0]]['name'] = ann[1]
    ind[ann[0]]['type'] = ann[2]
index.close()

head = next(table)
sys.stdout.write('name\ttype\t' + head)
for line in table:
    info = line.split('\t')
    if float(info[3]) > 0:
        sys.stdout.write(ind[info[1]]['name'] + '\t ' + ind[info[1]]['type'] + '\t' + line)
    else:
        sys.stderr.write('Skipped ' + ind[info[1]]['name'] + ' ' + ind[info[1]]['type'] + ' ' + info[1] + ' no reads!\n')
table.close()
