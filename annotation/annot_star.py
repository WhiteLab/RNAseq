#!/usr/bin/python
"""
Usage: annot_star.py <ct_file> <gtf>

Arguments:
<ct_file> read count file from htseq or STAR
<gtf>     gtf reference file

Options:
-h
"""
from docopt import docopt

args = docopt(__doc__)

import sys
import re

htseq = open(args['<ct_file>'], 'r')
hdict = {}
for line in htseq:
    line = line.rstrip('\n')
    data = line.split('\t')
    #    if data[0][0:2]=='DQ':
    #        pdb.set_trace()
    hdict[data[0]] = {}
    hdict[data[0]]['ct'] = '\t'.join(data[1:])
    hdict[data[0]]['f'] = 0
htseq.close()
gtf = open(args['<gtf>'], 'r')
sys.stdout.write('name\tid\ttype\tunstranded\tfirst_strand\tsecond_strand\n')
for line in gtf:
    if line[0] == '#':
        continue
    line = line.rstrip('\n')
    cols = line.split('\t')
    try:
        m = re.search('gene_id "(\S+)"; transcript_id "(\S+)";.*transcript_type "(\S+)";.*transcript_name "(\S+)";',
                      cols[-1])
        (gid, tid, ty, tn) = (m.group(1), m.group(2), m.group(3), m.group(4))
    #        if gid[0:2]=='DQ':
    #            pdb.set_trace()
    except:
        sys.stderr.write('Regex failed, skipping!\n')
        continue
    if gid in hdict and hdict[gid]['f'] == 0:
        sys.stdout.write(tn + '\t' + gid + '\t' + ty + '\t' + hdict[gid]['ct'] + '\n')
        hdict[gid]['f'] = 1
    else:
        if tid in hdict:
            sys.stdout.write(tn + '\t' + tid + '\t' + ty + '\t' + hdict[tid]['ct'] + '\n')
            hdict[tid]['f'] = 1
gtf.close()
for ids in hdict:
    if hdict[ids]['f'] == 0:
        sys.stdout.write('NA\t' + ids + '\tNA\t' + hdict[ids]['ct'] + '\n')
