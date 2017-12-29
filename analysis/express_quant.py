#!/usr/bin/env python3
import sys
import os
from utility.date_time import date_time
from utility.log import log
import subprocess
import json


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    return config_data['params']['stranded'], config_data['params']['strand'], config_data['tools']['express'],\
           config_data['refs']['express']


def express_quant(sample, config_file, ref_mnt, x, s):
    loc = sample + '.express.log'
    if os.path.isdir('LOGS'):
        loc = 'LOGS/' + loc
    (stranded, strand, express, transcriptome) = parse_config(config_file)

    transcriptome = ref_mnt + '/' + transcriptome
    if stranded == 'N':
        express_cmd = express + ' ' + transcriptome + ' ' + sample + '.merged.transcriptome.bam --no-update-check -m '\
                      + x + ' -s ' + s + ' --logtostderr 2>> ' + loc
    else:
        express_cmd = express + ' ' + transcriptome + ' ' + sample + '.merged.transcriptome.bam --no-update-check --'\
                      + strand + ' -m ' + x + ' -s ' + s + ' --logtostderr 2>> ' + loc
    log(loc, date_time() + express_cmd + '\n')
    check = subprocess.call(express_cmd, shell=True)

    rename_express_out = 'mv results.xprs ' + sample + '.express_quantification.txt; mv params.xprs ' + sample\
                         + '.params.xprs'
    check += subprocess.call(rename_express_out, shell=True)
    log(loc, date_time() + 'Completed qc.  Renaming files\n')
    return check

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')
    parser.add_argument('-x', '--mean', action='store', dest='x',
                        help='Mean insert size')
    parser.add_argument('-s', '--stdev', action='store', dest='s',
                        help='Standard deviation of insert size')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    sample = inputs.sample
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt
    x = inputs.x
    s = inputs.s
    express_quant(sample, config_file, ref_mnt, x, s)
