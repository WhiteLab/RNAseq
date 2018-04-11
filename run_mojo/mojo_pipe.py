#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
from utility.date_time import date_time
from utility.log import log
import subprocess
import json
from utility.set_acls import set_acls


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'], \
               config_data['tools']['mojo'], config_data['refs']['mojo_config'], config_data['params']['threads'], \
               config_data['params']['ram']
    except:
        try:
            sys.stderr.write(date_time() + 'Accessing keys failed.  Attempting to output current keys:\n')
            for key in config_data:
                sys.stderr.write(key + '\n')
                for subkey in config_data[key]:
                    sys.stderr.write(key + ":" + subkey + ":" + config_data[key][subkey] + '\n')
            exit(1)
        except:
            sys.stderr.write(date_time() + 'Could not read config file ' + json_config + '\n')
            exit(1)


def mojo_pipe(sample, config_file, fq1, fq2):
    (project_dir, project, align_dir, mojo, m_config, cores, mem) = parse_config(config_file)
    fq_dir = project + project + '/' + align_dir + '/' + sample + '/TRIMMED_FQ/'
    out_dir = project + project + '/' + align_dir + '/' + sample + '/MOJO/'
    os.chdir(fq_dir)
    os.mkdir(out_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run Chai\'s MOJO software on a single sample')
    parser.add_argument('-s', '--sample', action='store', dest='sample', help='Sample name')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-a', '--fastq1', action='store', dest='fq1',
                        help='Comma seprated list of read 1 fastq files')
    parser.add_argument('-b', '--fastq2', action='store', dest='fq2',
                        help='Comma seprated list of read 2 fastq files')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    mojo_pipe(inputs.sample, inputs.config_file, inputs.fq1, inputs.fq2)
