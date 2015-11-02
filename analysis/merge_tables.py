#!/usr/bin/env python
"""
Merges tables based on a single column, column denoted array-style. Designed for variant output - will merge gene,
position, and base change as name for row

Usage: ./merge_tables.py (<list> <col> <hflag> <tflag>) [options]

Arguments:
<list>   list of tables
<col>    index of column to merge, array style
<hflag>  1 if tables have header, 0 if not
<tflag>  1 if showing on-target only

Options:
-h
-s SUFFIX  if given, will omit from column name.  otherwise file name used as column header
-a ACCEPT  if given, a list of acceptable effects - i.e. NON_SYNONYMOUS
"""
import os
import sys

from docopt import docopt

args = docopt(__doc__)
tlist = open(args['<list>'])
col = int(args['<col>'])
hflag = int(args['<hflag>'])
tflag = int(args['<tflag>'])

sflag = 1
suffix = ''
if '-s' in args:
    suffix = args['-s']
else:
    sys.stderr.write('No suffix given.  Using file names as headers\n')
    sflag = 0
alist = {}
aflag = 0
if '-a' in args:
    aflag = 1
    list_h = open(args['-a'], 'r')
    for line in list_h:
        line = line.rstrip('\n')
        alist[line] = 1
new_tbl = {}
flist = []
# chr1	2488217	GTTxTGA	C	A	149	0	0.00%	99	5	4.81%	10000.00	NA	TNFRSF14	INTRON	CODING			283	KEEP	ON
temp = {}
vlist = []
for tbl in tlist:
    sys.stderr.write('Processing table ' + tbl)
    tbl = tbl.rstrip('\n')
    fh = open(tbl, 'r')
    tbl = os.path.basename(tbl)
    if sflag:
        tbl = tbl.replace(suffix, '')
    flist.append(tbl)
    if hflag:
        head = next(fh)
    new_tbl[tbl] = new_tbl[tbl] = {}
    for line in fh:
        line = line.rstrip('\n')
        if len(line) < 1:
            continue
        data = line.split('\t')
        if tflag and data[-1] == 'OFF':
            continue
        if aflag and data[14] not in alist:
            continue
        # pdb.set_trace()
        var = data[3] + '-' + data[4]
        row = '_'.join([data[13], data[0], data[1], var])
        if row not in temp:
            temp[row] = 1
            vlist.append(row)
        new_tbl[tbl][row] = data[col]
    fh.close()
tlist.close()
sys.stdout.write('Sample/variant' + '\t' + '\t'.join(flist) + '\n')
flist.sort()
vlist.sort()
for variant in vlist:
    sys.stdout.write(variant)
    for samp in flist:
        if variant in new_tbl[samp]:
            sys.stdout.write('\t' + new_tbl[samp][variant])
        else:
            sys.stdout.write('\t0')
    sys.stdout.write('\n')
