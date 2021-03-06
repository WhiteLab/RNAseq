#!/usr/bin/python
import sys
import os
import re
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
from utility.date_time import date_time
from utility.log import log


def align_stats(sample):
    # casual logging - look for a LOGS directory, otherwise assume current dir
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.aln.log'
    log(loc, date_time() + "Converting to table summary format\n")
    fh = open(sample + '/' + 'align_summary.txt', 'r')
    fo = open(sample + '.align.txt', 'w')
    fo.write(
        'Sample\tMean insert size estimate(10k reads)\tStd dev read insert size estimate(10 k reads)\tStarting left reads\t% mapped\tmultimapped(mm)\tgt 20 mm\tStarting right reads\t% mapped\t% mm\tgt 20 mm\tOverall map rate\tAligned pairs\t% mm\t% discordant\t% condordant\n' + sample + '\t')
    fi = open(sample + '_subset.insert_metrics.hist')
    for i in range(0, 7, 1):
        skip = next(fi)
    stats = next(fi)
    fi.close()
    stat = stats.split('\t')
    fo.write('\t'.join([str(int(float(stat[4]))), str(int(float(stat[5])))]))
    next(fh)
    lstart = next(fh)
    m = re.search('(\d+)\n$', lstart)
    fo.write('\t' + m.group(1))
    pct = next(fh)
    m = re.search('\(\s*(\S+) of input\)\n', pct)
    fo.write('\t' + m.group(1))
    mm = next(fh)
    m = re.search('\(\s*(\S+)\).*\((\d+) have >20\)\n', mm)
    fo.write('\t' + m.group(1) + '\t' + m.group(2))

    next(fh)
    rstart = next(fh)
    m = re.search('(\d+)\n$', rstart)
    fo.write('\t' + m.group(1))
    pct = next(fh)
    m = re.search('\(\s*(\S+) of input\)\n', pct)
    fo.write('\t' + m.group(1))
    mm = next(fh)
    m = re.search('\(\s*(\S+)\).*\((\d+) have >20\)\n', mm)
    fo.write('\t' + m.group(1) + '\t' + m.group(2))
    ovr = next(fh)
    m = re.search('\s*(^\S+)', ovr)
    fo.write('\t' + m.group(1))
    next(fh)

    aln = next(fh)
    m = re.search('(\d+)\n$', aln)
    fo.write('\t' + m.group(1))
    mm = next(fh)
    m = re.search('\(\s*(\S+)\) have', mm)
    fo.write('\t' + m.group(1))
    dc = next(fh)
    m = re.search('\(\s*(\S+)\) are', dc)
    fo.write('\t' + m.group(1))
    cc = next(fh)
    m = re.search('^\s*(\S+)', cc)
    fo.write('\t' + m.group(1) + '\n')
    fo.close
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Alignment summary report.  Converts tophat alignment summary to a table format')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/location name prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample) = (inputs.sample)
    align_stats(sample)
