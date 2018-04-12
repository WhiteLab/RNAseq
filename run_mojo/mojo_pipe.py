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
               config_data['params']['ram'], config_data['params']['user'], config_data['params']['group']
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
    (project_dir, project, align_dir, mojo, m_config, cores, mem, user, group) = parse_config(config_file)
    fq_dir = project_dir + project + '/' + align_dir + '/' + sample + '/TRIMMED_FQ/'
    out_dir = project_dir + project + '/' + align_dir + '/' + sample + '/MOJO/'
    os.mkdir(out_dir)
    loc = out_dir + sample + '.mojo_run.log'
    log(loc, date_time() + 'Made output directory ' + out_dir + '\n')
    log(loc, date_time() + 'Changing to fastq directory ' + fq_dir + '\n')
    os.chdir(fq_dir)
    run_mojo = mojo + ' --config ' + m_config + ' --sample_name ' + sample + ' --output_dir ' + out_dir + ' --fq1 ' \
               + fq1 + ' --fq2 ' + fq2 + ' --cores ' + cores + ' --mem ' + mem
    log(loc, date_time() + 'Running MOJO with command ' + run_mojo + '\n')
    try:
        subprocess.call(run_mojo, shell=True)
        log(loc, date_time() + 'MOJO complete! Setting acls\n')
        check = set_acls(out_dir, user, group)
        if check == 0:
            log(loc, date_time() + 'Setting acls complete.  Pipeline complete!\n')
        else:
            log(loc, date_time() + 'Setting acls failed.  Check logs!\n')
    except:
        sys.stderr.write(date_time() + 'MOJO failed!  Check logs in ' + loc + '\n')
        return 1


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
