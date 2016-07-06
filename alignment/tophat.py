#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log
import subprocess


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def tophat(tophat_tool, tx, bwt2_ref, end1, end2, x, s, sample, log_dir, th):
    loc = log_dir + sample + ".tophat.log"
    # fix to determine which phred score to use using HGAC date assigned 150409 and greater phred33, else hpred 64
    meta = sample.split('_')
    epoch = 150409
    tophat_cmd = tophat_tool + " --no-coverage-search --mate-inner-dist " + x + " --mate-std-dev " + s\
                 + " --num-threads " + th + " --library-type fr-firststrand --transcriptome-index " + tx + " -o " \
                 + sample + " " + bwt2_ref + " " + end1 + " " + end2 + " 2>> " + loc

    if len(meta) >= 2 and RepresentsInt(meta[1]) == True and int(meta[1]) < epoch:
        tophat_cmd = tophat_tool + " --no-coverage-search  --phred64-quals --mate-inner-dist " + x + " --mate-std-dev "\
                     + s + " --num-threads " + th + " --library-type fr-firststrand --transcriptome-index " + tx\
                     + " -o " + sample + " " + bwt2_ref + " " + end1 + " " + end2 + " 2>> " + loc
    log(loc, date_time() + tophat_cmd + "\n")
    call(tophat_cmd, shell=True)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='tophat paired-end alignment and transcript assembly module.')
    parser.add_argument('-t', '--tophat', action='store', dest='tophat_tool', help='Location of tophat alignment tool.')
    parser.add_argument('-tx', '--transcriptome', action='store', dest='tx', help='Location of pre-built transcriptome')
    parser.add_argument('-b', '--bwt2_reference', action='store', dest='bwt2_ref',
                        help='Location of bowtie2 reference file')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-x', '--mean', action='store', dest='x', help='Mean insert size')
    parser.add_argument('-sd', '--standard_deviation', action='store', dest='s',
                        help='Standard deviation of insert size')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-th', '--threads', action='store', dest='th', help='Number of threads')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (tophat_tool, tx, bwt2_ref, end1, end2, x, s, sample, log_dir, th) = (
    inputs.tophat_tool, inputs.tx, inputs.bwt2_ref, inputs.end1, inputs.end2, inputs.x, inputs.s, inputs.sample,
    inputs.log_dir, inputs.th)
    tophat(tophat_tool, tx, bwt2_ref, end1, end2, x, s, sample, log_dir, th)
