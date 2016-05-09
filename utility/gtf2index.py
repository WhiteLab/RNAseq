#!/usr/bin/env python

import sys
import re

gtf = sys.argv[1]
an_type = sys.argv[2]

fh = open(gtf, 'r')
print 'ID\tName\tType'

idf = an_type + '_id'
namef = an_type + '_name'
typef = an_type + '_type'
(id_name, name, type_name) = ('', '', '')
for line in fh:
    if line[0] != '#':
        info = line.rstrip('\n').split('\t')
        if info[2] == an_type:
            id_name = re.search(idf + ' "(\S+)"', info[8])
            name = re.search(namef + ' "(\S+)"', info[8])
            type = re.search(typef + ' "(\S+)"', info[8])
            print '\t'.join((id_name.group(1), name.group(1), type.group(1)))
fh.close()
