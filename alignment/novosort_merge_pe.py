#!/usr/bin/python
import json
import os
import re
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from subprocess import call
import subprocess
from picard_mark_dups import picard_mark_dups
from utility.job_manager import job_manager


def parse_config(config_file):
    config_data = json.load(open(config_file, 'r'))
    return config_data['tools']['novosort'], config_data['refs']['cont'], config_data['refs']['obj'], \
           config_data['params']['threads'], config_data['params']['ram'], config_data['params']['mark_dup']


def list_bam(cont, obj, sample, th, in_suffix):
    ct = 0
    list_cmd = '. /home/ubuntu/.novarc;swift list ' + cont + ' --prefix ' + obj + '/' + sample
    sys.stderr.write(date_time() + list_cmd + '\nGetting BAM list\n')
    flist = subprocess.check_output(list_cmd, shell=True)
    # Use to check on download status
    p = []
    bam_list = []
    for fn in re.findall('(.*)\n', flist):
        if re.match('^\S+_\d+' + in_suffix + '$', fn):
            sys.stderr.write(date_time() + 'Downloading relevant BAM file ' + fn + '\n')
            dl_cmd = '. /home/ubuntu/.novarc;swift download ' + cont + ' --skip-identical ' + fn + ' >> LOGS/' \
                     + sample + '.novosort_merge.log'
            p.append(dl_cmd)
            if fn[-3:] == 'bam':
                bam_list.append(fn)
                ct = ct + 1

    f = job_manager(p, th)
    if f == 0:
        sys.stderr.write(date_time() + 'BAM download complete for sample ' + sample + '\n')
        return bam_list, ct
    else:
        sys.stderr.write(date_time() + 'BAM download failed for sample ' + sample + '\n')
        exit(1)


def novosort_merge_pe(config_file, sample_list, in_suffix, out_suffix, sort_type):
    (novosort, cont, obj, th, ram, mflag) = parse_config(config_file)
    # gives some flexibility if giving a list of samples ot just a single one
    if os.path.isfile(sample_list):
        fh = open(sample_list, 'r')
    else:
        fh = []
        fh.append(sample_list)
    # create temp dir for sorting
    mk_temp = 'mkdir nova_temp'
    call(mk_temp, shell=True)
    for sample in fh:
        sample = sample.rstrip('\n')
        (bam_list, n) = list_bam(cont, obj, sample, th, in_suffix)
        if mflag == 'Y':
            sys.stderr.write(date_time() + 'Mark duplicates flag given.  Running picard on bams for sample' + sample
                             + '\n')
            for i in xrange(0, len(bam_list), 1):
                cur = bam_list[i].replace(in_suffix, '')
                check = picard_mark_dups(config_file, cur, 'LOGS/', in_suffix)
                if check != 0:
                    sys.stderr.write('Picard mark dup failed for ' + sample + '\n')
                    exit(1)

                bam_list[i] = bam_list[i].replace(in_suffix, '.dup_marked.bam')
        bam_string = " ".join(bam_list)
        final_bam = sample + out_suffix
        if sort_type == 'name':
            novosort_merge_pe_cmd = novosort + " -c " + th + " -m " + ram + "G  -o " + final_bam + ' -n -t' \
                                ' nova_temp ' + bam_string + ' 2>> LOGS/' + sample + '.novosort_merge.log'
        else:
            novosort_merge_pe_cmd = novosort + " -c " + th + " -m " + ram + "G  -o " + final_bam + ' -i -t' \
                                ' nova_tempz` ' + bam_string + ' 2>> LOGS/' + sample + '.novosort_merge.log'
        sys.stderr.write(date_time() + novosort_merge_pe_cmd + "\n")
        try:
            subprocess.check_output(novosort_merge_pe_cmd, shell=True)
        except:
            sys.stderr.write(date_time() + 'novosort failed for sample ' + sample + '\n')
            exit(1)
    rm_temp = 'rm -rf nova_temp'
    call(rm_temp, shell=True)
    sys.stderr.write(date_time() + 'Merge process complete\n')
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to merge BAM files module.')
    parser.add_argument('-sl', '--sample_list', action='store', dest='sample_list', help='Sample/project prefix list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and ref locations')
    parser.add_argument('-os', '--out_bam_suffix', action='store', dest='out_suffix', help='Suffix of output bam')
    parser.add_argument('-is', '--in_bam_suffix', action='store', dest='in_suffix', help='Suffix of input bam')
    parser.add_argument('-t', '--sort_type', action='store', dest='sort_type', help='Name or coordinate sort')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample_list, config_file, in_suffix, out_suffix, sort_type) = (inputs.sample_list, inputs.config_file,
                                                                inputs.in_suffix, inputs.out_suffix, inputs.sort_type)
    novosort_merge_pe(sample_list, config_file, in_suffix, out_suffix, sort_type)
