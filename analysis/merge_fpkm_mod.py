#!/usr/bin/python

import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time

flist = open(sys.argv[1], 'r')
data = {}
bids = []
for fn in flist:
    fn = fn.rstrip('\n')
    sys.stderr.write(date_time() + 'Processing file ' + fn + '\n')
    parts = fn.split('/')
    bid = parts[1]
    bids.append(bid)
    cur = open(fn, 'r')
    head = next(cur)
    for line in cur:
        line = line.rstrip('\n')
        datum = line.split('\t')
        # will only bother outputting transcripts with values > 0
        tx = datum[4]
        val = datum[9]
        if tx == '-':
            tx = datum[6]
        if float(val) > 0:
            if tx not in data:
                data[tx] = {}
            data[tx][bid] = val
    sys.stderr.write(date_time() + 'Completed processing file ' + fn + '\n')
    cur.close()
flist.close()
sys.stderr.write(date_time() + 'Outputting new master table\n')
sys.stdout.write('transcript/sample')
for bid in bids:
    sys.stdout.write('\t' + bid)
sys.stdout.write('\n')
for tx in data:
    sys.stdout.write(tx)
    for bid in bids:
        sys.stdout.write('\t')
        if bid in data[tx]:
            sys.stdout.write(data[tx][bid])
        else:
            sys.stdout.write('0')
    sys.stdout.write('\n')
sys.stderr.write(date_time() + 'Fin\n')
