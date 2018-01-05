#!/usr/bin/env python3

import json
import sys
import os
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
from utility.date_time import date_time
from utility.log import log
import subprocess


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['novosort'], config_data['tools']['java'], config_data['tools']['picard'], \
           config_data['refs']['project'], config_data['refs']['project_dir'], config_data['refs']['align_dir'], \
           config_data['params']['threads'], config_data['params']['ram'], \
           config_data['tools']['novo_merge_rmdup_slurm']


def list_bam(project, align, sample):
    bam_dir = '/cephfs/PROJECTS/' + project + '/' + align + '/' + sample + '/BAM/'
    find_bam_cmd = 'find ' + bam_dir + ' -name *.Aligned.toTranscriptome.out.bam'
    sys.stderr.write(date_time() + find_bam_cmd + '\nGetting BAM list\n')
    try:
        bam_find = subprocess.check_output(find_bam_cmd, shell=True).decode().rstrip('\n')
        bam_list = bam_find.split('\n')
        ct = len(bam_list)

        return bam_list, ct
    except:
        sys.stderr.write(date_time() + 'No bams found for ' + sample + '\n')
        exit(1)


def novosort_merge_pe(config_file, sample_list):
    fh = open(sample_list, 'r')
    (novosort, java_tool, picard_tool, project, project_dir, align, threads, ram, novo_merge_rmdup_slurm) \
        = parse_config(config_file)

    for sample in fh:
        sample = sample.rstrip('\n')
        loc = '../LOGS/' + sample + '.novosort_merge.log'
        job_loc = sample + '.novosort_merge.log'
        (bam_list, n) = list_bam(project, align, sample)
        bam_string = " ".join(bam_list)
        cur_dir = project_dir + project + '/' + align + '/' + sample + '/BAM/'
        os.chdir(cur_dir)
        out_bam = sample + '.merged.transcriptome.bam'
        if n > 1:
            batch = 'sbatch -c ' + threads + ' --mem ' + ram + ' -o ' + job_loc + ' --export=novosort="' \
                    + novosort + '",threads="' + threads + '",ram="' + ram + 'G",out_bam="' + out_bam \
                    + '",bam_string="' + bam_string + '",loc="' + loc + '"' + ' ' + novo_merge_rmdup_slurm
            log(loc, date_time() + 'Submitting merge bam job for sample ' + batch + "\n")
            subprocess.call(batch, shell=True)


        else:

                link_bam = 'ln -s ' + bam_list[0] + ' ' + sample + '.merged.transcriptome.bam;'
                log(loc, date_time() + 'Creating symlink for merged final bam since only one exists\n'
                    + link_bam + '\n')
                subprocess.call(link_bam, shell=True)

    sys.stderr.write(date_time() + 'Merged file request submitted and processed, check logs.\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample/project prefix list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file) = (inputs.sample_list, inputs.config_file)
    novosort_merge_pe(config_file, sample_list)
