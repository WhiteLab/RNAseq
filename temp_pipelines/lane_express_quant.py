#!/usr/bin/env python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import os
import json
from analysis.express_quant import parse_config
from utility.job_manager import job_manager


def lane_express_quant(bams, config_file, ref_mnt):
    (stranded, strand, express, transcriptome) = parse_config(config_file)
    transcriptome = ref_mnt + '/' + transcriptome
    job_list = []
    th = '4'
    for bam in open(bams):
        bam = bam.rstrip('\n')
        (bam_dir, root) = (os.path.dirname(bam), os.path.basename(bam.replace('.Aligned.toTranscriptome.out.bam', '')))
        parts = root.split('_')
        qc_file = 'ALIGN_RNASEQ/' + parts[0] + '/QC/' + root + '.qc_stats.json'
        qc_data = json.loads(open(qc_file, 'r').read())
        (x, s) = (str(int(round(qc_data['picard_stats']['x_ins_size']))),
                  str(int(round(qc_data['picard_stats']['s_ins_size']))))
        wd = bam_dir + '/' + root + '/'
        loc = wd + root + '.log'
        express_cmd = 'mkdir ' + wd + ';'
        if stranded == 'N':
            express_cmd += express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + bam_dir + ' -m '\
                          + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'
        else:
            express_cmd += express + ' ' + transcriptome + ' ' + bam + ' --no-update-check -o ' + bam_dir + ' --'\
                          + strand + ' -m ' + x + ' -s ' + s + ' --logtostderr 2>> ' + loc + ';'

        express_cmd += 'mv ' + wd + 'results.xprs ' + wd + root + '.express_quantification.txt; mv ' + wd \
                       + 'params.xprs ' + wd + root + '.params.xprs;'
        job_list.append(express_cmd)
        job_manager(job_list, th)
        return 0

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Quantify transcripts using STAR output bam')
    parser.add_argument('-b', '--bams', action='store', dest='bams', help='bam file list')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    bams = inputs.bams
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt
    lane_express_quant(bams, config_file, ref_mnt)
