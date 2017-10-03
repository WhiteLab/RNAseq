#!/usr/bin/env python
"""
Combines express output tables, collapsing on gene symbol by effective count

Usage: collapse_merge_express_by_gene.py <ct_list> <type_list> <field> <round>

Arguments:
  <ct_list> list of express count files
  <type_list> list of transcript types to accept. Type None not to use one
  <field> name of field to collapse on or default for est_counts
  <round> \'y\' to round values after collapsing

Options:
  -h
"""
from docopt import docopt

args = docopt(__doc__)
import sys
import os
import re

s_list = []
g_dict = {}
g_list = []
tbl_dict = {}
ty_flag = 'N'
ty_dict = {}

if args['<type_list>'] != 'None':
    ty_flag = 'Y'
    for line in open(args['<type_list>']):
        line = line.rstrip('\n')
        ty_dict[line] = 1

for report in open(args['<ct_list>']):
    report = report.rstrip('\n')
    fn = os.path.basename(report)
    parts = fn.split('.')
    s_list.append(parts[0])
    fh = open(report)
    head = next(fh)
    head_vals = head.rstrip('\n').split('\t')
    f = 9
    if args['<field>'] != 'default':
        flag = 0
        for i in xrange(len(head_vals)):
            if head_vals[i] == args['<field>']:
                f = i
                flag = 1
                break
        if flag == 0:
            sys.stderr.write('Field name ' + args['<field>'] + ' not found.  Using default!\n')

    for line in fh:
        info = line.rstrip('\n').split()
        gene = re.sub('-\d+$', '', info[0])
        if ty_flag == 'N' or info[1] in ty_dict:
            if gene not in g_dict:
                g_list.append(gene)
                g_dict[gene] = 1
                tbl_dict[gene] = {}
            if parts[0] not in tbl_dict[gene]:
                tbl_dict[gene][parts[0]] = 0.0
            tbl_dict[gene][parts[0]] += float(info[f])
    fh.close()

print 'Sample/Gene\t' + '\t'.join(s_list)
if args['<round>'] == 'y':
    for gene in g_list:
        sys.stdout.write(gene)
        for samp in s_list:
            if samp in tbl_dict[gene]:
                sys.stdout.write('\t' + str(int(round(tbl_dict[gene][samp]))))
            else:
                sys.stdout.write('\t0')
        print
else:
    for gene in g_list:
        sys.stdout.write(gene)
        for samp in s_list:
            if samp in tbl_dict[gene]:
                sys.stdout.write('\t' + str(tbl_dict[gene][samp]))
            else:
                sys.stdout.write('\t0')
        print
