#!/usr/bin/env python
import sys
import os
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
sys.path.append('/home/ubuntu/TOOLS/Scripts/annotation')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from log import log
from job_manager import job_manager
from novosort_merge_pe import novosort_merge_pe
import subprocess
import json


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['refs']['cont'], config_data['refs']['obj'], \
               config_data['refs']['capture_flag'], config_data['params']['capture_intvl'], \
               config_data['refs']['cap_bed'], config_data['tools']['bedtools'], \
               config_data['tools']['samtools'], config_data['tools']['java'], config_data['tools']['gatk'], \
               config_data['params']['threads'], config_data['refs']['samtools']
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


def filter_bam(bedtools, cap_bed, sample_list, out_suffix, th):
    job_list = []
    for sample in sample_list:
        bam = sample + out_suffix
        filt_cmd = bedtools + ' intersect -abam ' + bam + ' -b ' + cap_bed + ' -wa -header > ' + bam + 'FILTERED; mv '\
                   + bam + ' FILTERED ' + bam + ';'
        job_list.append(filt_cmd)
    check = job_manager(job_list, th)
    if check == 0:
        return 0
    else:
        return 1


def splitNtrim(java, gatk, sample_list, out_suffix, fasta, th):
    cmd_list = []
    for sample in sample_list:
        bam = sample_list + out_suffix
        split_cmd = java + ' -jar ' + gatk + ' -T SplitNCigarReads -R ' + fasta + ' -I ' + bam \
                    + ' -o ' + sample + '.merged.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60' \
                                        ' -U ALLOW_N_CIGAR_READS'
        cmd_list.append(split_cmd)
    rflag = job_manager(cmd_list, th)
    if rflag == 0:
        return 0
    else:
        return 1


def base_recal(java, gatk, sample_list, th, fasta):
    for sample in sample_list:
        bam = sample + '.merged.split.bam'
        recal_cmd = java + ' -jar ' + gatk + ' -nct ' + th + ' -T BaseRecalibrator -I ' + bam + ' -o ' + sample \
                    + '_recal_data.table -R ' + fasta
        rflag = subprocess.call(recal_cmd, shell=True)
        if rflag != 0:
            return 1
        new_bam_cmd = java + ' -jar ' + gatk + ' -nct ' + th + ' -T PrintReads -I ' + bam + ' -BQSR ' + sample \
                      + '_recal_data.table -o ' + sample + '.recalibrated.bam'
        rflag = subprocess.call(new_bam_cmd, shell=True)
        if rflag == 0:
            rm_old_bam = 'rm ' + bam
            subprocess.call(rm_old_bam, shell=True)
            return 0
        else:
            return 1


def gatk_call(sample_pairs, config_file, ref_mnt):
    mk_dir = 'mkdir BAM LOGS ANALYSIS ANNOTATION'
    subprocess.call(mk_dir, shell=True)
    (cont, obj, cflag, cintvl, cap_bed, bedtools, samtools, java, gatk, th, fasta) = parse_config(config_file)
    cintvl = ref_mnt + '/' + cintvl
    cap_bed = ref_mnt + '/' + cap_bed
    fasta = ref_mnt + '/' + fasta

    sample_list = 'sample_list.txt'
    slist = []
    fh = open(sample_pairs, 'r')
    sl = open(sample_list, 'w')
    temp = {}
    for line in fh:
        cur = line.rstrip('\n').split('\t')
        if cur[1] not in temp:
            sl.write(cur[1] + '\n')
            temp[cur[1]] = 1
            slist.append(cur[1])
        if cur[2] not in temp:
            sl.write(cur[2] + '\n')
            temp[cur[2]] = 1
            slist.append(cur[2])
    sl.close()
    fh .close()
    del temp
    in_suffix = '.Aligned.sortedByCoord.out.bam'
    out_suffix = '.merged.sortedByCoord.bam'
    sort_type = 'coordinate'
    check = novosort_merge_pe(sample_list, config_file, in_suffix, out_suffix, sort_type)
    if check == 0:
        sys.stderr.write(date_time() + 'File download and merge complete!\n')
        # rm unmerged bams, no longer needed
        rm_bam = 'rm -rf ' + obj
        subprocess.call(rm_bam, shell=True)
    else:
        sys.stderr.write(date_time() + 'File download and merge failed!\n')
        exit(1)

    if cflag == 'Y':
        sys.stderr.write(date_time() + 'Capture flag detected.  Filtering bams by capture regions\n')
        rflag = filter_bam(bedtools, cap_bed, slist, out_suffix, th)
        if rflag != 0:
            sys.stderr.write(date_time() + 'Filter bam failed!\n')
            exit(1)
        else:
            sys.stderr.write(date_time() + 'Filtering complete\n')
    check = splitNtrim(java, gatk, slist, out_suffix, fasta, th)
    if check != 0:
        sys.stderr.write(date_time() + 'Split n trim failed\n')
        exit(1)

    check = base_recal(java, gatk, slist, th, fasta)
    if check != 0:
        sys.stderr.write(date_time() + 'Base recal failed\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Use GATK haplotype caller best practices for mutation calls')
    parser.add_argument('-sp', '--pairs', action='store', dest='pairs', help='Sample piars')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    pairs = inputs.pairs
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt

    gatk_call(pairs, config_file, ref_mnt)