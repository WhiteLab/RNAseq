#!/usr/bin/python

"""
Collapses fastq file to fasta file with read count as part of sequence header

Usage: fq2collapsedfa.py <fastq> 

Arguments:
<fastq> fastq file

Options:
-h
"""
from docopt import docopt

args = docopt(__doc__)
import sys
import gzip

if args['<fastq>'][-2:] == 'gz':
    fq = gzip.open(args['<fastq>'])
else:
    fq = open(args['<fastq>'])
seq_dict = {}
x = 1
m = 1000000
for line in fq:
    seq = next(fq)
    if (x % m) == 0:
        sys.stderr.write('Processing seq ' + str(x) + '\n')
    seq = seq.rstrip('\n')
    if seq not in seq_dict:
        seq_dict[seq] = 0
    seq_dict[seq] += 1
    skip = next(fq)
    skip = next(fq)
    x += 1
fq.close()
i = 1
for seq in sorted(seq_dict, key=seq_dict.get, reverse=True):
    sys.stdout.write('>seq' + str(i) + '|' + str(seq_dict[seq]) + '\n' + seq + '\n')
    i += 1
