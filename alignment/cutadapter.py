#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
import json
from utility.date_time import date_time
from subprocess import call
from utility.log import log


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['cutadapt'], config_data['params']['threads'], config_data['params']['minlen'], \
           config_data['params']['r1adapt'], config_data['params']['r2adapt'], config_data['params']['r1trim'], \
           config_data['params']['r2trim'], config_data['params']['qual'], config_data['params']['mqual']


def cutadapter(sample, end1, end2, config_file):
    # casual logging - look for a LOGS directory, otherwise assume current dir
    log_dir = './'
    # designed to be run in a subdirectory, keep original file names
    sf1 = end1
    sf2 = end2
    end1 = os.path.basename(sf1)
    end2 = os.path.basename(sf2)
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.cutadapt.log'
    (cutadapt_tool, threads, minlen, r1adapt, r2adapt, r1trim, r2trim, qual, mqual) = parse_config(config_file)
    cut_th = threads
    if int(cut_th) >= 4:
        cut_th = str(int(threads) - 2)

    cutadapt_cmd = cutadapt_tool + ' -j ' + cut_th + ' -m ' + minlen + ' --quality-base=' + qual + ' -q ' + mqual \
                   + ' -a ' + r1adapt + ' -A ' + r2adapt + ' -u ' + r1trim + ' -U ' + r2trim + ' -o ' + end1 \
                   + ' -p ' + end2 + ' ' + sf1 + ' ' + sf2 + ' >> ' + loc + ' 2>> ' + loc
    log(loc, date_time() + cutadapt_cmd + "\n")
    call(cutadapt_cmd, shell=True)
    return 0


if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description='cutadapt module.  Removes 3\' adapters and trims bases if necessary.Also can enforce minimum read length - 20 recommended for tophat')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/location name prefix')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample, end1, end2, config_file) = (inputs.sample, inputs.end1, inputs.end2, inputs.config_file)
    cutadapter(sample, end1, end2, config_file)
