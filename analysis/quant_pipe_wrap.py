#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
from utility.date_time import date_time
from subprocess import call
import json
from statistics import mean


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'], \
               config_data['tools']['quant_slurm_wrap'], config_data['tools']['quant_pipe'], \
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


def quant_pipe_wrap(lane, config_file):
    (project_dir, project, align_dir, quant_slurm_wrap, quant_pipe, cores, mem) = parse_config(config_file)

    fh = open(lane, 'r')
    for sample in fh:
        (bnid, stype, lanes) = sample.rstrip('\n').split('\t')
        sys.stderr.write(date_time() + 'Getting ready to quantify ' + bnid + '\n')
        lanes = lanes.split(', ')
        x = []
        s = []
        f = 0
        for cur_lane in lanes:
            qc_file = project_dir + project + '/' + align_dir + '/' + bnid + '/QC/' + bnid + '_' + cur_lane \
                      + '.qc_stats.json'
            if os.path.isfile(qc_file):
                qc_data = json.loads(open(qc_file, 'r').read())
                x.append(float(qc_data['picard_stats']['x_ins_size']))
                s.append(float(qc_data['picard_stats']['s_ins_size']))
            else:
                sys.stderr.write(date_time() + 'Could not find qc file ' + qc_file + '!\nCheck inputs, skipping '
                                                                                     'sample ' + bnid + '\n')
                f = 1
        if f == 0:
            cur_mean = int(round(mean(x)))
            cur_std = int(round(mean(s)))
            job_name = 'rnaseq-quant_' + bnid
            job_log = bnid + '.quant.log'
            batch = 'sbatch -J ' + job_name + ' -c ' + cores + ' --mem ' + mem + ' -o ' + job_log \
                    + ' --export=quant="' + quant_pipe + '",bnid="' + bnid + '",j="' + config_file \
                    + '",x="' + str(cur_mean) + '",s="' + str(cur_std) + '" ' + quant_slurm_wrap
            sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
            try:
                call(batch, shell=True)
            except:
                sys.stderr.write(date_time() + 'Batch submission for ' + bnid + ' failed! Check logs!\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-l', '--lane', action='store', dest='lane', help='Lane list file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    lane = inputs.lane
    config_file = inputs.config_file

    quant_pipe_wrap(lane, config_file)
