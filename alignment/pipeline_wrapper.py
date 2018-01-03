#!/usr/bin/env python3
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
import re
import argparse
import json
from utility.date_time import date_time
from utility.find_project_files import find_project_files
from subprocess import call
from utility.log import log

parser = argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially.')
parser.add_argument('-f', '--file', action='store', dest='fn',
                    help='File with bionimbus ID, seqtype and sample lane list')
parser.add_argument('-j', '--json', action='store', dest='config_file',
                    help='JSON config file with tools, references, and data storage locations')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

inputs = parser.parse_args()
fh = open(inputs.fn, 'r')


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    # skip pdx flag allows to pick up from post-filtering step to save substantial time of just human align need repeat
    return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'], \
           config_data['refs']['config'], config_data['params']['threads'], config_data['params']['ram'], \
           config_data['tools']['align_pipe'], config_data['tools']['slurm_wrap']


(project_dir, project, align_dir, pipe_cfg, cores, mem, align_pipe, slurm_wrap) = parse_config(inputs.config_file)
cwd = project_dir + project
if not os.path.isdir(cwd):
    sys.stderr.write(date_time() + 'Could not find working directory ' + cwd
                     + '. Ensure correct project was set in config\n')
    exit(1)
for line in fh:
    line = line.rstrip('\n')
    (bnid, seqtype, lane_csv) = line.split('\t')
    cur_dir = cwd + '/RAW/' + bnid
    for lane in lane_csv.split(', '):
        loc = bnid + '_' + lane + '.log'
        file_prefix = bnid + '_' + lane
        (contents, seqfile, sf1, sf2) = ('', [], '', '')
        # attempt to find sequencing files
        try:
            sys.stderr.write(date_time() + 'Searching for sequencing files related to bnid ' + bnid + ' in lane '
                             + lane + '\n')
            contents = find_project_files(cur_dir, file_prefix)
            # standard file naming should work with this
            seqfile = re.findall('(\S+[sequence|f*q]*\.gz)', contents)
            # need to sort since output may not be alphanumerically sorted
            seqfile.sort()
            sf1 = seqfile[0]
            sf2 = seqfile[1]
        except:
            sys.stderr.write(date_time() + 'Getting sequencing files for ' + bnid + ' lane ' + lane
                             + ' failed.  Moving on\n')
            continue
        log(loc, date_time() + 'Running pipeline process for lane ' + lane + '\n')
        # check class status flag
        batch = 'sbatch -c ' + cores + ' --mem ' + mem + ' -o ' + loc \
                + ' --export=pipeline="' + align_pipe + '",f1="' + sf1 + '",f2="' + sf2 + '",j="' + pipe_cfg + '"' \
                + ' ' + slurm_wrap
        sys.stderr.write(date_time() + 'Submitting job ' + batch + '\n')
        try:
            call(batch, shell=True)
        except:
            sys.stderr.write(date_time() + 'Batch submission for ' + lane + ' failed! Check logs!\n')

sys.stderr.write(date_time() + "Jobs submitted.  Check logs for any errors\n")