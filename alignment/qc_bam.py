#!/usr/bin/env python3

import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
import os
from utility.date_time import date_time
import subprocess
import json
from utility.log import log


def parse_config(json_config):
    config_data = json.loads(open(json_config, 'r').read())
    return config_data['tools']['java'], config_data['params']['ram'], config_data['tools']['picard'],\
           config_data['refs']['refFlat'], config_data['refs']['rRNA_intervals'], config_data['params']['strand'], \
           config_data['params']['threads']


def qc_bam(sample, config_file, ref_mnt):
    # job_list = []
    loc = sample + '.bam_qc.log'
    if os.path.isdir('LOGS'):
        loc = 'LOGS/' + loc
    (java, ram, picard, refFlat, intervals, strand, threads) = parse_config(config_file)
    # recalc ram to be a bit lower
    ram = str(int(round(int(ram) * 0.75)))

    st_dict = {'N': 'NONE', 'fr-stranded': 'FIRST_READ_TRANSCRIPTION_STRAND',
               'rf-stranded': 'SECOND_READ_TRANSCRIPTION_STRAND'}

    refFlat = ref_mnt + '/' + refFlat
    intervals = ref_mnt + '/' + intervals
    picard_cmd = java + ' -Xmx' + ram + 'g -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=' + threads +  \
                 ' -XX:MaxGCPauseMillis=10000 -jar ' + picard + ' CollectRnaSeqMetrics REF_FLAT=' + refFlat \
                 + ' STRAND=' + st_dict[strand] + ' CHART=' + sample + '.pos_v_cov.pdf I=' + sample \
                 + '.Aligned.sortedByCoord.out.bam O=' + sample + '.picard_RNAseq_qc.txt RIBOSOMAL_INTERVALS=' \
                 + intervals + ' VALIDATION_STRINGENCY=SILENT 2>> ' + loc + ' >> ' + loc
    # job_list.append(picard_cmd)
    log(loc, date_time() + picard_cmd + '\n')
    subprocess.call(picard_cmd, shell=True)
    return 0


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
    qc_bam(sample, config_file, ref_mnt)


if __name__ == "__main__":
    main()
