#!/usr/bin/env python

# helper script for converting and sorting a table

import sys
import re
# row or column
sort_type = sys.argv[1]
# for now, csv or tsv
space_type = sys.argv[2]
# where to start sorting, rest will be fixed!
start = int(sys.argv[3])
table = sys.argv[4]


def make_csv(data):
    data = re.sub('\t', ',', data)
    return data


def make_tsv(data):
    data = re.sub(',', '\t', data)
    return data

fh = open(table)
new_table = {}
if sort_type == 'row':
    if space_type == 'csv':
        for i in xrange(start):
            line = next(fh)
            line = make_csv(line)
            sys.stdout.write(line)
    elif space_type == 'tsv':
        for i in xrange(start):
            line = next(fh)
            line = make_tsv(line)
            sys.stdout.write(line)
    else:
        for i in xrange(start):
            line = next(fh)
            sys.stdout.writable(line)
    if space_type == 'csv':
        for line in fh:
            line = line.rstrip('\n')
            line = make_csv(line)
            new_table[line] = 1
    elif space_type == 'tsv':
        for line in fh:
            line = line.rstrip('\n')
            line = make_csv(line)
            new_table[line] = 1
    else:
        for line in fh:
            line = line.rstrip('\n')
            new_table[line] = 1

    fh.close()
    for line in sorted(new_table):
        sys.stdout.write(line + '\n')

if sort_type == 'col':
    head = next(fh)
    sort_key = []
    if space_type == 'csv':
        head = head.rstrip('\n').split('\t')
        unsort_head = []
        to_sort_head = head[start:]
        for i in xrange(start):
            unsort_head.append(head[i])
        sort_key = sorted(range(len(to_sort_head)), key=lambda k: to_sort_head[k])
        sys.stdout.write(','.join(unsort_head))
        for i in sort_key:
            sys.stdout.write(',' + to_sort_head[i])
        print
        for line in fh:
            head = line.rstrip('\n').split('\t')
            unsort_head = []
            to_sort_head = head[start:]
            for i in xrange(start):
                unsort_head.append(head[i])
            sys.stdout.write(','.join(unsort_head))
            for i in sort_key:
                sys.stdout.write(',' + to_sort_head[i])
            print
    fh.close()


