#!/usr/bin/env python
"""
Collapses tables based on a gene name - can choose to make values binary

Usage: ./collapse_var.py (<table> <bflag>)

Arguments:
<table>    table
<bflag>    1 to make values binary, 0 to sum

Options:
-h
"""
import re
import sys

from docopt import docopt

args = docopt(__doc__)
fh = open(args['<table>'], 'r')
bflag = args['<bflag>']

head = next(fh)
head = head.rstrip('\n')
hlist = head.split('\t')
print head
genes = {}

for line in fh:
    line = line.rstrip('\n')
    data = line.split('\t')
    m = re.search('(\S+)_chr.*', data[0])
    try:
        gene = m.group(1)
    except:
        sys.stderr.write('Failed to find gene while processing ' + line + '\n')
        exit(1)
    if gene not in genes:
        genes[gene] = {}
    for i in xrange(1, len(data), 1):
        if bflag == '1' and int(data[i]) > 0:
            genes[gene][hlist[i]] = 1
        else:
            if hlist[i] not in genes[gene]:
                genes[gene][hlist[i]] = 0
            genes[gene][hlist[i]] += int(data[i])
fh.close()
for gene in sorted(genes.keys()):
    sys.stdout.write(gene)
    for i in xrange(1, len(hlist), 1):
        if hlist[i] in genes[gene]:
            sys.stdout.write('\t' + str(genes[gene][hlist[i]]))
        else:
            sys.stdout.write('\t0')
    sys.stdout.write('\n')
