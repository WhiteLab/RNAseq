#!/usr/bin/python
import json
import os
import re
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
import subprocess
from job_manager import job_manager


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['novosort'], config_data['refs']['cont'], config_data['refs']['obj'],\
           config_data['params']['threads'], config_data['params']['ram'], config_data['tools']['samtools']


def list_bam(cont, obj, sample, wait, th):
    ct = 0
    list_cmd = '. /home/ubuntu/.novarc;swift list ' + cont + ' --prefix ' + obj + '/' + sample
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    p = []
    bam_list = []
    bai_list = []
    for fn in re.findall('(.*)\n', flist):
        if re.match('^\S+_\d+\.Aligned.toTranscriptome.out.bam$', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + fn
            p.append(dl_cmd)
            if fn[-3:] == 'bam':
                bam_list.append(fn)
                ct = ct + 1
            #else:
            #    bai_list.append(fn)
    f = job_manager(p, th)
    if f == 0:
        sys.stderr.write(date_time() + 'BAM download complete\n')
        return bam_list, bai_list, ct
    else:
        sys.stderr.write(date_time() + 'BAM download failed\n')
        exit(1)


def novosort_merge_pe(config_file, sample_list, wait):
    (novosort, cont, obj, th, ram, samtools) = parse_config(config_file)
    # gives some flexibility if giving a list of samples ot just a single one
    if os.path.isfile(sample_list):
        fh = open(sample_list, 'r')
    else:
        fh = []
        fh.append(sample_list)
    for sample in fh:
        sample = sample.rstrip('\n')
        (bam_list, bai_list, n) = list_bam(cont, obj, sample, wait, th)
        bam_string = ",".join(bam_list)
        final_bam = sample + '.merged.final.bam'
        if n > 1:
            novosort_merge_pe_cmd = novosort + " --c " + th + " --m " + ram + "G --assumesorted --output " + final_bam \
                                    + ' --index --md --tmpdir ./TMP ' + bam_string
            sys.stderr.write(date_time() + novosort_merge_pe_cmd + "\n")
            try:
                subprocess.check_output(novosort_merge_pe_cmd, shell=True)
            except:
                sys.stderr.write(date_time() + 'novosort failed for sample ' + sample + '\n')
                exit(1)
        else:
            rename_bam = 'cp ' + bam_list[0] + ' ' + sample + '.merged.final.bam; ' + samtools + ' index ' + final_bam
            sys.stderr.write(date_time() + rename_bam + ' Only one associated bam file, renaming\n')
            call(rename_bam, shell=True)
    sys.stderr.write(date_time() + 'Merge process complete\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample/project prefix list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-w', '--wait', action='store', dest='wait',
                        help='Wait time to download bam files.  900 (seconds) recommended')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file, wait) = (inputs.sample_list, inputs.config_file, inputs.wait)
    novosort_merge_pe(config_file, sample_list, wait)
