#!/usr/bin/env python

import sys


cat_file = open(sys.argv[1])
cont_file = open(sys.argv[2])
table = open(sys.argv[3])

cat_dict = {}
cont_dict = {}
cat_list = []
cont_list = []

for line in cat_file:
    line = line .rstrip('\n')
    cat_dict[line] = 0
    cat_list.append(line)
cat_file.close()

for line in cont_file:
    line = line .rstrip('\n')
    cont_dict[line] = 0
    cont_list.append(line)
cont_file.close()

head = next(table)
head = head.rstrip('\n').split('\t')
new_head = []
for i in xrange(len(head)):
    if head[i] in cat_dict:
        cat_dict[head[i]] = i
        new_head.append(head[i])
    elif head[i] in cont_dict:
        cont_dict[head[i]] = i
        new_head.append(head[i])

print '\t'.join(new_head)
meta_dict = {}
for line in table:
    info = line.rstrip('\n').split('\t')
    if info[0] not in meta_dict:
        meta_dict[info[0]] = {}
        meta_dict[info[0]]['cat'] = []
        for cat in cat_list:
            meta_dict[info[0]]['cat'].append(info[cat_dict[cat]])
        meta_dict[info[0]]['cont'] = {}
        meta_dict[info[0]]['cont'][cont_list[0]] = 0
        for i in xrange(1, len(cont_list), 1):
            meta_dict[info[0]]['cont'][cont_list[i]] = [0, 0]
    meta_dict[info[0]]['cont'][cont_list[0]] += int(info[cont_dict[cont_list[0]]].replace(',', ''))
    for i in xrange(1, len(cont_list), 1):
        value = info[cont_dict[cont_list[i]]].replace('%', '')
        meta_dict[info[0]]['cont'][cont_list[i]][0] += float(value)
        meta_dict[info[0]]['cont'][cont_list[i]][1] += 1
table.close()

for bnid in meta_dict:
    sys.stdout.write('\t'.join(meta_dict[bnid]['cat']) + '\t'
                     + str(meta_dict[bnid]['cont'][cont_list[0]]))
    for i in xrange(1, len(cont_list), 1):
        avg = (meta_dict[bnid]['cont'][cont_list[i]][0]/meta_dict[bnid]['cont'][cont_list[i]][1])
        sys.stdout.write('\t' + str(avg))
    print


