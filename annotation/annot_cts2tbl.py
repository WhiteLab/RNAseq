#!/usr/bin/python
"""
Converts read count outputs to table
Usage: annot_cts2table.py <ct_file_list> <gtf>

Arguments:
<ct_file_list> read count file list from htseq or STAR
<gtf>     gtf reference file

Options:
-h
"""
from docopt import docopt
args = docopt(__doc__)

import sys
import re
import pdb
import os

ct_list=open(args['<ct_file_list>'],'r')
ctdict={}
slist = []
for line in ct_list:
    line = line.rstrip('\n')
    RG = os.path.basename(line)
    samp = RG.split('_')
    sn = samp[0]
    slist.append(sn)
    cur = open(line,'r')
    for cts in cur:
        cts = cts.rstrip('\n')
        data = cts.split('\t')
        #    if data[0][0:2]=='DQ':
        #        pdb.set_trace()
        if data[0] not in ctdict:
            ctdict[data[0]]={}
        ctdict[data[0]][sn] = {}
        ctdict[data[0]][sn]['ct']='\t'.join(data[1:])
        ctdict[data[0]]['f']=0
    cur.close()
ct_list.close()
gtf=open(args['<gtf>'],'r')
#sys.stdout.write('name\tid\ttype\tunstranded\tfirst_strand\tsecond_strand\n')
cout = ('unstranded','first_strand','second_strand')
fh = {}
for strand in cout:
    fh[strand] = open(strand + '_cts.tsv','w')
    fh[strand].write('name\tid\ttype')
    for samp in slist:
        fh[strand].write('\t' + samp)
for strand in cout:
    fh[strand].write('\n')
for line in gtf:
    if line[0] == '#':
        continue
    line = line.rstrip('\n')
    cols = line.split('\t')
    try:
        m=re.search('gene_id "(\S+)"; transcript_id "(\S+)";.*transcript_type "(\S+)";.*transcript_name "(\S+)";',cols[-1])
        (gid,tid,ty,tn)=(m.group(1),m.group(2),m.group(3),m.group(4))
#        if gid[0:2]=='DQ':
#            pdb.set_trace()
    except:
        sys.stderr.write('Regex failed, skipping!\n')
        continue
    if gid in ctdict and ctdict[gid]['f']==0:
        for i in xrange(0,len(cout),1):
            fh[cout[i]].write(tn + '\t' + gid + '\t' + ty)
            for samp in slist:
                vals = ctdict[gid][samp]['ct'].split('\t')
                fh[cout[i]].write('\t' + vals[i])
        for i in xrange(0,len(cout),1):
            fh[cout[i]].write('\n')
        ctdict[gid]['f']=1
    else:
        if tid in ctdict:
            fh[cout[i]].write(tn + '\t' + gid + '\t' + ty)
            for samp in slist:
                vals = ctdict[tid][samp]['ct'].split('\t')
                fh[cout[i]].write('\t' + vals[i])
            for i in xrange(0,len(cout),1):
                fh[cout[i]].write('\n')
            ctdict[tid]['f']=1
gtf.close()
for ids in ctdict:
    if ctdict[ids]['f']==0:
        #sys.stdout.write('NA\t' + ids + '\tNA\t' + ctdict[ids]['ct'] + '\n') 
        for i in xrange(0,len(cout),1):
            fh[cout[i]].write('NA\t' + ids + '\tNA')
            for samp in slist:
                vals = ctdict[ids][samp]['ct'].split('\t')
                fh[cout[i]].write('\t' + vals[i])
        for i in xrange(0,len(cout),1):
            fh[cout[i]].write('\n')
