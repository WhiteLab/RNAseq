#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/RNAseq')
from utility.date_time import date_time
from temp_pipelines.downsample_bam import downsample_bam
import subprocess
import os
import json


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['samtools'], config_data['params']['threads']


def get_from_depth(qc_root, depth):
    qc_data = json.loads(open(qc_root + '.qc_stats.json', 'r').read())
    (start_rp, passed_rp) = (qc_data['cutadapt_stats']['starting_read_pairs'], qc_data['cutadapt_stats']['rp_pass'])
    start_rp = int(start_rp.replace(',', ''))
    passed_rp = float(passed_rp.replace('%', ''))/100
    passed_rd_ct = start_rp * passed_rp
    frac = round(int(depth)/passed_rd_ct, 2)
    return str(frac)


def downsample_pipe(bam_list, config_file, depth):
    (samtools, threads) = parse_config(config_file)
    for bam in open(bam_list):
        sys.stderr.write(date_time() + 'Setting up for ' + bam)
        bam = bam.rstrip('\n')
        bam_root = bam.replace('.bam', '')
        bam_dir = os.path.dirname(bam)
        # cleanup sample name, sloppy i know!
        qc_root = bam_root.replace('BAMS', 'QC', 1)
        qc_root = qc_root.replace('.srt', '', 1)
        qc_root = qc_root.replace('.Aligned.toTranscriptome.out', '', 1)
        sys.stderr.write(date_time() + 'Calculating downsample fraction\n')
        frac = get_from_depth(qc_root, depth)
        # submit to job queue
        downsample_bam(samtools, bam, frac, bam_dir, threads)
        #cmd = ' '.join(('sbatch', '-c', threads, '--oversubscribe', downsample_bam, '-b ', bam, '-f', frac, '-o ',
        #                bam_dir, '-t', threads, '-s', samtools))
        sys.stderr.write(date_time() + 'Submitting to queue ' + bam + '\n')
        #subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Downsample bams to desired read depth')
    parser.add_argument('-l', '--bam_list', action='store', dest='bam_list', help='List of bam files')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-d', '--depth', action='store', dest='depth',
                        help='desired read depth')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    bam_list = inputs.bam_list
    config_file = inputs.config_file
    depth = inputs.depth
    downsample_pipe(bam_list, config_file, depth)
