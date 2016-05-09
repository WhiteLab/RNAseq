#!/usr/bin/env python
import sys

cfile = open(sys.argv[1], 'r')
genome = open(sys.argv[2], 'r')

clist = []

for line in cfile:
    line = line.rstrip('\n')
    clist.append(line)
cfile.close()

gdict = {}
cur = ''
for line in genome:
    if line[0] == ">":
        cur = line[1:]
        gdict[cur] = []
    else:
        gdict[cur].append(line)
genome.close()

for chrom in clist:
    sys.stdout.wite('>' + chrom + '\n' + ''.join(gdict[chrom]))
