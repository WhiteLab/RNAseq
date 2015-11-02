#!/usr/bin/python

import sys

fh = open(sys.argv[1], 'r')
cur = []
for line in fh:
    line = line.rstrip('\n')
    cur.append(line)
rev = list(reversed(cur))
fh.close()
print '\n'.join(rev)
