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

sys.stdout.write('name\ttype\t')
head = next(table)
for line in table:
    info = line.split('\t')
    sys.stdout.write(ind[info[0]]['name'] + '\t ' + ind[info[0]]['type'] + '\t' + line)

table.close()