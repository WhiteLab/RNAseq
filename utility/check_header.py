#!/usr/bin/python
import sys
import gzip
import pdb

f1=gzip.open(sys.argv[1],'rb')
f2=gzip.open(sys.argv[2],'rb')

m = 4000000

i=1
c=10000

while i < m:
    if (i % c) ==0:
        sys.stderr.write('At line ' + str(i) + '\n')
    line1=next(f1)
    line2=next(f2)
    id1=line1.split()
    id2=line2.split()
    seq1=next(f1)
    seq2=next(f2)
#    pdb.set_trace()
    if id1[0] != id2[0]:
        sys.stderr.write(line1 + line2 + 'are discordant at ' + str(i) + ' with sequences\n' + seq1 + seq2)
        exit(1)
    skip1=next(f1)
    skip1=next(f1)
    skip1=next(f2)
    skip1=next(f2)
    i+=4
f1.close()
f2.close()
    
