#!/usr/bin/env python3

import sys
import os
import json
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from subprocess import call
from utility.log import log


def parse_config(config_file):
    config_data = json.load(open(config_file, 'r'))
    return config_data['tools']['java'], config_data['tools']['picard'], config_data['params']['ram']


def picard_mark_dups(config_file, sample, log_dir, suffix):
    root = os.path.basename(sample)
    loc = log_dir + root + ".picard.mark_dup.log"
    (java_tool, picard_tool, mem) = parse_config(config_file)
    picard_tmp = 'picard_tmp'
    picard_mark_dups_cmd = 'mkdir ' + picard_tmp + ';' + java_tool + " -Djava.io.tmpdir=" + picard_tmp + " -Xmx" \
                           + mem + "g -jar " + picard_tool + " MarkDuplicates I=" + sample + suffix + " O=" + sample \
                           + ".dup_marked.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=" + sample \
                           + ".output.metrics > " + loc + " 2>&1; rm -rf " + picard_tmp
    log(loc, date_time() + picard_mark_dups_cmd + "\n")
    check = call(picard_mark_dups_cmd, shell=True)
    if check == 0:
        return 0
    else:
        return 1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Picard mark duplicates module.  Needed for gatk haplotype caller preprocessing.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-s', '--suffix', action='store', dest='suffix', help='Bam file suffix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample, log_dir, suffix) = (inputs.config_file, inputs.sample,
                                                         inputs.log_dir, inputs.suffix)
    picard_mark_dups(config_file, sample, log_dir, suffix)
