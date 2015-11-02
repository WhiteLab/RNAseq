#!/usr/bin/python
import sys

fh = open(sys.argv[1], 'r')
ind = {}
for line in fh:
    line = line.rstrip('\n')
    cols = line.split('\t')
    ind[cols[0]] = '\t'.join(cols[1:])
fh.close()
kout = open(sys.argv[2])
head = next(kout)
head = head.rstrip('\n')
sys.stdout.write(head + '\tref_id\ttx_name\ttype\n')
for line in kout:
    line = line.rstrip('\n')
    cols = line.split('\t')
    sys.stdout.write(line + '\t' + ind[cols[0]] + '\n')
kout.close()
