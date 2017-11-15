#!/usr/bin/env python
import sys
sys.path.append('/cephfs/users/mbrown/RNAseq')
import os
import json
from analysis.express_quant import parse_config
from subprocess import call


def lane_express_quant(bams, config_file):
    (stranded, strand, express, transcriptome) = parse_config(config_file)
    for bam in open(bams):
        bam = bam.rstrip('\n')
        (bam_dir, root) = (os.path.dirname(bam), os.path.basename(bam.replace('.Aligned.toTranscriptome.out.*', '')))
        qc_dir = bam_dir.replace('BAMS', 'QC')
        qc_file = qc_dir + '/' + root + '.qc_stats.json'
        qc_data = json.loads(open(qc_file, 'r').read())
        (x, s) = (str(int(round(float(qc_data['picard_stats']['x_ins_size'])))),
                  str(int(round(float(qc_data['picard_stats']['s_ins_size'])))))
        wd = qc_dir + '/' + root + '/'
        loc = wd + root + '.log'
        express_cmd = 'mkdir ' + wd + ';'
        if stranded == 'N':
            express_cmd += express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + wd + ' -m '\
                          + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'
        else:
            express_cmd += express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + wd + ' --'\
                          + strand + ' -m ' + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'

        express_cmd += 'mv ' + wd + 'results.xprs ' + wd + root + '.express_quantification.txt; mv ' + wd \
                       + 'params.xprs ' + wd + root + '.params.xprs;'
        call('sbatch -c 4 ' + express_cmd, shell=True)

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
