#!/usr/bin/env python

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
from utility.date_time import date_time
from utility.log import log
from express_quant import express_quant
from alignment.novosort_merge_pe import novosort_merge_pe
from annotation.annot_express import annot_express
import subprocess
import json
from statistics import mean


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['refs']['cont'], config_data['refs']['obj'], config_data['refs']['tx_index'], \
               config_data['params']['capture_flag'], config_data['refs']['cap_targets'],\
               config_data['tools']['samtools'], config_data['tools']['filter_bam']
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


def upload_special(bnid, cont, obj):
    src_cmd = '. ~/.novarc;'
    bam = obj + '/' + bnid + '/BAMS/' + bnid + '.merged.transcriptome.bam'

    ONE_GB = 1073741824
    up_bam = src_cmd + ' swift upload -S ' + str(ONE_GB) + ' ' + cont + ' ' + bam
    check = subprocess.call(up_bam, shell=True)
    if check != 0:
        sys.stderr.write('Could not upload bam file. Command given: ' + up_bam + '\n')
    report_dir = obj + '/' + bnid + '/REPORTS/'
    up_reports = src_cmd + ' swift upload ' + cont + ' ' + report_dir
    check = subprocess.call(up_reports, shell=True)
    if check != 0:
        sys.stderr.write('Could not upload report file. Command given: ' + up_reports + '\n')
    log_dir = obj + '/' + bnid + '/LOGS/'
    up_logs = src_cmd + ' swift upload ' + cont + ' ' + log_dir
    check = subprocess.call(up_logs, shell=True)
    if check != 0:
        sys.stderr.write('Could not upload report file. Command given: ' + up_logs + '\n')


def quant_pipe(lane, config_file, ref_mnt):
    src_cmd = '. ~/.novarc;'
    (cont, obj, tx_index, cflag, ctargets, samtools, filter_bam) = parse_config(config_file)
    mk_tmp = 'mkdir nova_temp'
    subprocess.call(mk_tmp, shell=True)
    if not os.path.isdir('BAMS'):
        mk_bam_dir = 'mkdir BAMS'
        subprocess.call(mk_bam_dir, shell=True)
    if not os.path.isdir('REPORTS'):
        mk_rpt_dir = 'mkdir REPORTS'
        subprocess.call(mk_rpt_dir, shell=True)
    fh = open(lane, 'r')
    for sample in fh:
        (bnid, stype, lanes) = sample.rstrip('\n').split('\t')
        log_dir = 'LOGS/'
        if not os.path.isdir(log_dir):
            mk_log_dir = 'mkdir ' + log_dir
            subprocess.call(mk_log_dir, shell=True)
        loc = log_dir + bnid + '.quantification_pipe.log'
        lanes = lanes.split(', ')
        x = []
        s = []
        for cur_lane in lanes:
            qc_file = obj + '/' + bnid + '/QC/' + bnid + '_' + cur_lane + '.qc_stats.json'
            get_qc = src_cmd + 'swift download ' + cont + ' ' + qc_file
            subprocess.call(get_qc, shell=True)
            qc_data = json.loads(open(qc_file, 'r').read())
            x.append(float(qc_data['picard_stats']['x_ins_size']))
            s.append(float(qc_data['picard_stats']['s_ins_size']))

        cur_mean = int(round(mean(x)))
        cur_std = int(round(mean(s)))
        in_suffix = '.Aligned.toTranscriptome.out.bam'
        out_suffix = '.merged.transcriptome.bam'
        sort_type = 'name'
        check = novosort_merge_pe(inputs.config_file, bnid, in_suffix, out_suffix, sort_type)
        # remove downloaded bam directory after merging
        rm_bam_dir = 'rm -rf ' + obj + '/' + bnid + '/BAMS'
        subprocess.call(rm_bam_dir, shell=True)
        if check != 0:
            log(loc, date_time() + 'Merge of RNAseq bams for ' + bnid + ' failed.  Please check logs\n')
            exit(1)
        if cflag == 'Y':
            mbam = sample + out_suffix
            out_bam = sample + '.merged.capture_filtered.transcriptome.bam'
            cmd = samtools + ' view -h ' + mbam + ' | ' + filter_bam + ' -i ' + ctargets + ' | ' + samtools\
                  + ' view -bSh - > ' + out_bam + '; rm ' + mbam + '; mv ' + out_bam + ' ' + mbam
            log(loc, date_time() + 'Capture method flag is Y, pre-filtering transcript hits\n' + cmd + '\n')
        check = express_quant(bnid, inputs.config_file, ref_mnt, str(cur_mean), str(cur_std))
        if check != 0:
            log(loc, date_time() + 'Quantification of RNA failed.  Please check logs\n')
            exit(1)
        # annotate with gene symbol and type, ignore 0 values for ease of use
        check = annot_express(tx_index, bnid, ref_mnt)
        if check != 0:
            log(loc, date_time() + 'Annotation of eXpress file.  Please check logs\n')
            exit(1)
        mv_cmd = 'mv *.bam BAMS/; mkdir REPORTS; mv *xpr* REPORTS/; mv BAMS REPORTS LOGS ' + obj + '/' + bnid
        subprocess.call(mv_cmd, shell=True)
        sys.stderr.write(date_time() + 'Uploading merged bam and quant files for ' + bnid + '\n')
        upload_special(bnid, cont, obj)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-l', '--lane', action='store', dest='lane', help='Lane list file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    lane = inputs.lane
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt

    quant_pipe(lane, config_file, ref_mnt)
