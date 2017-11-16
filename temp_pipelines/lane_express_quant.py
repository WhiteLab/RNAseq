#!/usr/bin/env python
import sys
sys.path.append('/cephfs/users/mbrown/RNAseq')
import os
import json
import re
from subprocess import call
from utility.date_time import date_time


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    return config_data['params']['stranded'], config_data['params']['strand'], config_data['tools']['express'],\
           config_data['tools']['express_sl'], config_data['refs']['express']


def lane_express_quant(bams, config_file):
    (stranded, strand, express, express_sl, transcriptome) = parse_config(config_file)
    for bam in open(bams):
        bam = bam.rstrip('\n')
        bam_dir = os.path.dirname(bam)
        root = os.path.basename(re.sub('.Aligned.toTranscriptome.out.*', '', bam))
        qc_dir = bam_dir.replace('BAMS', 'QC')
        qc_file = qc_dir + '/' + root + '.qc_stats.json'
        qc_data = json.loads(open(qc_file, 'r').read())
        (x, s) = (str(int(round(float(qc_data['picard_stats']['x_ins_size'])))),
                  str(int(round(float(qc_data['picard_stats']['s_ins_size'])))))
        wd = qc_dir + '/' + root + '/'
        loc = wd + root + '.log'
        express_cmd = 'mkdir ' + wd + ';'
        call(express_cmd, shell=True)
        sys.stderr.write(date_time() + 'Created dir ' + wd + ' to quantify ' + bam + '\n' + express_cmd + '\n')
        if stranded == 'N':
            express_cmd = express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + wd + ' -m '\
                          + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'
        else:
            express_cmd = 'sbatch -c 4 --export=express="' + express + '",transcriptome="' + transcriptome + '",bam="' \
                          + bam + '",wd="' + wd + '",strand="' + strand + '",x="' + x + '",s="' + s + '",loc="' + loc \
                          + '",root="' + root + '" ' + express_sl
            # express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + wd + ' --'\
            #              + strand + ' -m ' + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'

            # express_cmd += 'mv ' + wd + 'results.xprs ' + wd + root + '.express_quantification.txt; mv ' + wd \
            #               + 'params.xprs ' + wd + root + '.params.xprs;'
        sys.stderr.write(date_time() + 'Submitting quantification job\n' + express_cmd + '\n')
        call(express_cmd, shell=True)

    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-b', '--bams', action='store', dest='bams', help='bam file list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    bams = inputs.bams
    config_file = inputs.config_file
    lane_express_quant(bams, config_file)
