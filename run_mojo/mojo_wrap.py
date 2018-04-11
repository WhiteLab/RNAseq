#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
from utility.date_time import date_time
from subprocess import call
import json


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['tools']['slurm_wrap'], config_data['tools']['mojo_pipe'], \
               config_data['params']['threads'], config_data['params']['ram']
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


def mojo_wrap(lane, config_file):
    (mojo_wrap, mojo_pipe, cores, mem) = parse_config(config_file)
    fh = open(lane, 'r')
    for sample in fh:
        (bnid, stype, lanes) = sample.rstrip('\n').split('\t')
        # fq_dir = project_dir + project + '/' + align_dir + '/' + bnid + '/TRIMMED_FQ/'
        fq1 = []
        fq2 = []
        for lane in lanes.split(', '):
            fq1.append(bnid + '_' + lane + '_1_sequence.txt.gz')
            fq2.append(bnid + '_' + lane + '_2_sequence.txt.gz')
        slurm_cmd = 'sbatch -J ' + bnid + '-MOJO -o ' + bnid + '_mojo.log -c ' + cores + ' --mem ' + mem \
                    + 'G export=mojo="' + mojo_pipe + '",j="' + config_file + '",a="' + ','.join(fq1) + '",b="' \
                    + ','.join(fq2) + '" ' + mojo_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + slurm_cmd + '\n')
        try:
            call(slurm_cmd, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + bnid + ' failed! Check logs!\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Search for fusion junctions using MOJO')
    parser.add_argument('-l', '--lane', action='store', dest='lane', help='Lane list file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    lane = inputs.lane
    config_file = inputs.config_file

    mojo_wrap(lane, config_file)
