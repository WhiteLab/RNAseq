#!/usr/bin/env python

import sys
sys.path.append('/cephfs/users/mbrown/RNAseq')
from utility.date_time import date_time
import subprocess
import os


def downsample_bam(samtools, bam, frac, out_dir, th):
    out_root = os.path.basename(bam.replace('.bam', ''))
    cmd = samtools + ' view --threads ' + th + ' -b ' + bam + ' -s ' + frac + ' > ' + out_dir + '/' + out_root \
                                          + '_subsample_' + frac + '.bam'
    sys.stderr.write(date_time() + 'Downsampling ' + bam + '\n' + cmd + '\n')
    subprocess.call(cmd, shell=True)
    sys.stderr.write(date_time() + 'process complete!\n')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Downsample bams for RNAseq anaylsis')
    parser.add_argument('-s', '--samtools', action='store', dest='samtools', help='samtools location')
    parser.add_argument('-b', '--bam_file', action='store', dest='bam', help='bam file name')
    parser.add_argument('-f', '--fraction', action='store', dest='frac',
                        help='Relative frequency of total reads to output')
    parser.add_argument('-o', '--output_dir', action='store', dest='out_dir',
                        help='Location to output resulting bam')
    parser.add_argument('-t', '--threads', action='store', dest='threads',
                        help='Num threads to use')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (samtools, bam, frac, out_dir, th) = (inputs.samtools, inputs.bam, inputs.frac, inputs.out_dir, inputs.threads)
    downsample_bam(samtools, bam, frac, out_dir, th)
