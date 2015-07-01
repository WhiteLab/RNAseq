#!/usr/bin/python
import sys
import re
import pdb

gtf=open(sys.argv[1])
ind={}
for line in gtf:
    if line[0] == '#':
        continue
    line = line.rstrip('\n')
    cols = line.split('\t')
    try:
        m=re.search('gene_biotype "(\S+)";.*transcript_id "(\S+)"; transcript_name "(\S+)"',cols[-1])
        (ty,tid,tn)=(m.group(1),m.group(2),m.group(3))
    except:
        m=re.search('transcript_id "(\S+)";.*transcript_type "(\S+)";.*transcript_name "(\S+)";',cols[-1])
        (tid,ty,tn)=(m.group(1),m.group(2),m.group(3))

    if tid not in ind:
        ind[tid]={}
        ind[tid]['t']=ty
        ind[tid]['tn']=tn
gtf.close()

fa=open(sys.argv[2])
for line in fa:
    if line[0] == '>':
        line=line.rstrip('\n')
        ids=line.split()
        try:
            eid=ids[1]
            sys.stdout.write(ids[0][1:] + '\t' + eid + '\t' + ind[eid]['tn'] + '\t' + ind[eid]['t'] + '\n')
        except:
            sys.stderr.write('Processing ' + line + ' threw an error\n')
            try:
                eid=line[1:]
                sys.stdout.write(eid + '\t' + eid + '\t' + ind[eid]['tn'] + '\t' + ind[eid]['t'] + '\n')
            except:
                pdb.set_trace()
                sys.stderr.write('Processing ' + line + ' threw an error\n')
                sys.stderr.write(eid + '\n')
                exit(1)
fa.close()
