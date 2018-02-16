#!/usr/bin/env python3
"""
Combines express output tables, collapsing on gene symbol by effective count

Usage: collapse_merge_express.py <ct_list> <type_list> <field> <round> <c_flag>

Arguments:
  <ct_list> list of express count files
  <type_list> list of transcript types to accept. Type None not to use one
  <field> name of field to collapse on or default for est_counts
  <round> \'y\' to round values after collapsing
  <c_flag> \'y\' to collapse by gene or just merge by transcript and ENSEMBL id

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


def by_gene(tbl_dict, g_dict, g_list, ty_flag, ty_dict, fh, f):
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
    return tbl_dict, g_dict, g_list


def by_tx(tbl_dict, g_dict, g_list, ty_flag, ty_dict, fh, f):
    for line in fh:
        info = line.rstrip('\n').split()
        gene = info[0] + '\t' + info[1]
        if ty_flag == 'N' or info[1] in ty_dict:
            if gene not in g_dict:
                g_list.append(gene)
                g_dict[gene] = 1
                tbl_dict[gene] = {}
            if parts[0] not in tbl_dict[gene]:
                tbl_dict[gene][parts[0]] = 0.0
            tbl_dict[gene][parts[0]] += float(info[f])
    fh.close()
    return tbl_dict, g_dict, g_list


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
        for i in range(len(head_vals)):
            if head_vals[i] == args['<field>']:
                f = i
                flag = 1
                break
        if flag == 0:
            sys.stderr.write('Field name ' + args['<field>'] + ' not found.  Using default!\n')
    if args['<c_flag>'] == 'y':
        (tbl_dict, g_dict, g_list) = by_gene(tbl_dict, g_dict, g_list, ty_flag, ty_dict, fh, f)
    else:
        (tbl_dict, g_dict, g_list) = by_tx(tbl_dict, g_dict, g_list, ty_flag, ty_dict, fh, f)

if args['<c_flag>'] == 'y':
    print ('Gene\t' + '\t'.join(s_list) + '\n')
else:
    print ('Tx_name\ttx_id\t' + '\t'.join(s_list) + '\n')
if args['<round>'] == 'y':
    for gene in g_list:
        sys.stdout.write(gene)
        for samp in s_list:
            if samp in tbl_dict[gene]:
                sys.stdout.write('\t' + str(int(round(tbl_dict[gene][samp]))))
            else:
                sys.stdout.write('\t0')
        print ('\n')
else:
    for gene in g_list:
        sys.stdout.write(gene)
        for samp in s_list:
            if samp in tbl_dict[gene]:
                sys.stdout.write('\t' + str(tbl_dict[gene][samp]))
            else:
                sys.stdout.write('\t0')
        print ('\n')
