#!/usr/bin/env python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
import json
from utility.job_manager import job_manager

def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['refs']['cont'], config_data['refs']['obj']


def qc_bam_pipe(sample, config_file, ref_mnt):
    (cont, obj) = parse_config(config_file)
    job_list = []
    log_dir = 'LOGS/'
    src_cmd = '. ~/.novarc;'
    for samp in open(sample):
        samp = samp.rstrip('\n')
        parts = samp.split('_')
        bam = samp + '.Aligned.sortedByCoord.out.bam'
        dl_list = (log_dir + samp + '.cutadapt.log', log_dir + samp + '.Log.final.out', 'QC/' + samp
                   + '_subset.insert_metrics.hist', 'QC/' + samp + '_1_sequence_fastqc/fastqc_data.txt', 'BAMS/'
                   + bam)
        dl_cmd = ''
        prefix = obj + '/' + parts[0] + '/'
        for fn in dl_list:
            dl_cmd += src_cmd + 'swift download ' + cont + ' ' + prefix + fn + ';'
        mv_cmd = 'mv ' + prefix + dl_list[2] + ' .;mv ' + prefix + dl_list[4] + ' .;'
        qc_cmd = '~/TOOLS/Scripts/alignment/qc_bam.py -sa ' + samp + ' -j ' + config_file + ' -m ' + ref_mnt + ';'
        rm_cmd = 'rm ' + bam + ';'
        full_cmd = dl_cmd + mv_cmd + qc_cmd + rm_cmd
        job_list.append(full_cmd)
    job_manager(job_list, 4)



def main():
    import argparse
    parser = argparse.ArgumentParser(description='Bam qc from STAR output')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    sample = inputs.sample
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt
    qc_bam_pipe(sample, config_file, ref_mnt)


if __name__ == "__main__":
    main()
