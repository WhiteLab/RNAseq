#!/usr/bin/env python3

import sys
import os
import json
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
from utility.date_time import date_time
from alignment.novosort_merge_pe import novosort_merge_pe


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['project'], config_data['refs']['project_dir'], config_data['refs']['align_dir']


def list_bam(project, project_dir, align_dir, sample):
    project_dir + project
    bam = project_dir + project + '/' + align_dir + '/' + sample + '/BAMS/' + sample + '.merged.transcriptome.bam'
    check_file = os.path.isfile(bam)
    if not check_file:
        sys.stderr.write(date_time() + 'Merged bam ' + bam + ' not found.\n')
    return check_file


def check_for_merged_bams(config_file, sample_list):
    fh = open(sample_list, 'r')
    (project, project_dir, align_dir) = parse_config(config_file)
    missing = []
    for sample in fh:
        sample = sample.rstrip('\n')
        check = list_bam(project, project_dir, align_dir, sample)
        if not check:
            missing.append(sample)
    return missing


def run_novosort(config_file, sample_list):
        check = novosort_merge_pe(config_file, sample_list)
        if check == 0:
            sys.stderr.write(date_time() + 'File merge complete!\n')

        else:
            sys.stderr.write(date_time() + 'File download and merge failed.\n')
            exit(1)


def preprocess_bams(config_file, lane_list):
    # create sample list
    sample_list = 'sample_list.txt'
    fh = open(lane_list, 'r')
    sl = open(sample_list, 'w')
    temp = {}
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        if cur[0] not in temp:
            sl.write(cur[0] + '\n')
            temp[cur[0]] = 1
    sl.close()
    fh .close()
    miss_list = check_for_merged_bams(config_file, sample_list)
    if len(miss_list) > 0:
        sys.stderr.write(date_time() + 'Missing files detected, merging lane files\n')
        temp_fn = 'temp_samp_list.txt'
        temp_fh = open(temp_fn, 'w')
        temp_fh.write('\n'.join(miss_list))
        temp_fh.close()
        run_novosort(config_file, temp_fn)
    else:
        sys.stderr.write(date_time() + 'All bams found. Ready for next step!\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Pre-variant calling step to merge all bam files for each sample '
                                                 'before running')
    parser.add_argument('-l', '--lane-list', action='store', dest='lane_list',
                        help='Tumor/normal sample pair list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (lane_list, config_file) = (inputs.lane_list, inputs.config_file)
    preprocess_bams(config_file, lane_list)
