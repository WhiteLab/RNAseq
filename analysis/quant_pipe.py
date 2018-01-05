#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
from utility.date_time import date_time
from utility.log import log
from express_quant import express_quant
from annotation.annot_express import annot_express
import subprocess
import json
from utility.set_acls import set_acls


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    try:
        return config_data['refs']['project_dir'], config_data['refs']['project'], config_data['refs']['align_dir'], \
               config_data['refs']['tx_index'], config_data['params']['capture_flag'], \
               config_data['refs']['cap_targets'], config_data['tools']['samtools'], \
               config_data['tools']['filter_bam'], config_data['params']['user'], config_file['params']['group']
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


def quant_pipe(sample, x, s, config_file):
    (project_dir, project, align_dir, tx_index, cflag, ctargets, samtools, filter_bam, user, group) \
        = parse_config(config_file)
    cwd = project_dir + project + '/' + align_dir + '/' + sample
    if os.path.isdir(cwd):
        sys.stderr.write(date_time() + 'Changing to working directory ' + cwd)
        os.chdir(cwd)
    else:
        sys.stderr.write(date_time() + 'Current working dir ' + cwd + ' Doesn\'t exist! Check inputs\n')
    bam_dir = 'BAMS/'
    report_dir = 'REPORTS/'
    if not os.path.isdir(report_dir):
        mk_rpt_dir = 'mkdir REPORTS'
        subprocess.call(mk_rpt_dir, shell=True)

    log_dir = 'LOGS/'
    if not os.path.isdir(log_dir):
        mk_log_dir = 'mkdir ' + log_dir
        subprocess.call(mk_log_dir, shell=True)
    loc = log_dir + sample + '.quantification_pipe.log'
    out_suffix = '.merged.transcriptome.bam'

    if cflag == 'Y':
        mbam = bam_dir + sample + out_suffix
        out_bam = bam_dir + sample + '.merged.capture_filtered.transcriptome.bam'
        cmd = samtools + ' view -h ' + mbam + ' | ' + filter_bam + ' -i ' + ctargets + ' | ' + samtools\
              + ' view -bSh - > ' + out_bam + '; rm ' + mbam + '; mv ' + out_bam + ' ' + mbam
        log(loc, date_time() + 'Capture method flag is Y, pre-filtering transcript hits\n' + cmd + '\n')
        check = subprocess.call(cmd, shell=True)
        if check != 0:
            sys.stderr.write('Prefilter hits failed for sample ' + sample + '\n')
            exit(1)
    check = express_quant(sample, inputs.config_file, str(x), str(s))
    if check != 0:
        log(loc, date_time() + 'Quantification of RNA failed.  Please check logs\n')
        exit(1)
    # annotate with gene symbol and type, ignore 0 values for ease of use
    check = annot_express(tx_index, sample)
    if check != 0:
        log(loc, date_time() + 'Annotation of eXpress file failed.  Please check logs\n')
        exit(1)
    mv_cmd = 'mv *xpr* REPORTS/;'
    log(loc, date_time() + 'Organizing outputs' + mv_cmd + '\n')
    subprocess.call(mv_cmd, shell=True)
    set_acls(report_dir, user, group)
    set_acls(bam_dir, user, group)
    set_acls(log_dir, user, group)
    sys.stderr.write('Quant pipe complete for ' + sample + ', check logs\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-b', '--sample', action='store', dest='sample', help='Sample ID')
    parser.add_argument('-x', '--mean', action='store', dest='x', help='Mean insert size rounded')
    parser.add_argument('-s', '--std', action='store', dest='s', help='Standard devation insert size rounded')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    sample = inputs.sample
    x = inputs.x
    s = inputs.s
    config_file = inputs.config_file

    quant_pipe(sample, x, s, config_file)
