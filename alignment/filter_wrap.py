#!/usr/bin/env python
import sys

from utility.date_time import date_time
from utility.log import log
import subprocess


def filter_wrap(mmu_filter, star_tool, genome_ref, end1, end2, sample, log_dir, threads, novosort):
    meta = sample.split('_')
    RGRP = "ID:" + sample + "\tLB:" + meta[0] + "\tPU:" + meta[4] + "\tSM:" + meta[0] + "\tPL:illumina"
    loc = log_dir + sample + ".mmu.star.pe.log"
    mk_srt_tmp = 'mkdir TMP'
    subprocess.call(mk_srt_tmp, shell=True)
    threads = str(int(threads) - 2)
    star_cmd = "(" + star_tool + " --runMode alignReads --outSAMattrRGline " + RGRP \
            + " --outFileNamePrefix " + sample + ".mmu_filt. --runThreadN " + threads + " --genomeDir " + genome_ref\
            + " --readFilesIn " + end1 + " " + end2 + " --readFilesCommand zcat --outSAMtype BAM Unsorted --outStd " \
            "BAM_Unsorted --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 " \
            "--alignSJDBoverhangMin 1 --outFilterMismatchNmax 0" + " --alignIntronMin 20 --alignIntronMax 1000000 " \
            "--alignMatesGapMax 1000000 --outSAMunmapped Within 2>> " + loc + "  | " + novosort + " - -n -c 2 -m 2G -t"\
            " TMP 2>> " + loc + " | tee " + sample + ".mmu.nsrt.bam | python " + mmu_filter + " -s " + \
            sample + " -n 0 -t RNA | gzip -4 -c - > " + sample + "_1.filtered.fq.gz;) 2>&1 | gzip -4 -c - > " + sample \
            + "_2.filtered.fq.gz"

    log(loc, date_time() + star_cmd + '\n')
    try:
        subprocess.call(star_cmd, shell=True)
    except:
        log(loc, date_time() + 'Star alignment and filter against against mouse genome failed\n')
        exit(1)
    log(loc, date_time() + 'Filtering completed, replacing fastq file\n')
    rn_fq = 'mv ' + sample + '_1.filtered.fq.gz ' + end1 + '; mv ' + sample + '_2.filtered.fq.gz ' + end2 \
            + ';rm -rf TMP'
    check = subprocess.call(rn_fq, shell=True)
    if check != 0:
        log(loc, date_time() + 'File rename failed\n' + rn_fq + '\n')
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mouse read filter alignment module. Fairly PDX-specific process.')
    parser.add_argument('-mmu', '--mmu_filter', action='store', dest='mmu_filter',
                        help='Location of bam filter tool.')
    parser.add_argument('-s', '--star', action='store', dest='star_tool',
                        help='Location of STAR alignment tool.  Version 2.4.2a preferred.')
    parser.add_argument('-g', '--genome_reference', action='store', dest='genome_ref', help='Location of star genome'
                                                                                            ' reference file')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-t', '--threads', action='store', dest='threads',
                        help='Number of threads to use.  8 recommended for standard vm')
    parser.add_argument('-n', '--novosort', action='store', dest='novosort',
                        help='Location of novosort to name sort read output')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (mmu_filter, star_tool, genome_ref, end1, end2, sample, log_dir, threads, novosort) = (
        inputs.mmu_filter, inputs.star_tool, inputs.genome_ref, inputs.end1, inputs.end2, inputs.sample,
        inputs.log_dir, inputs.threads, inputs.novosort)
    filter_wrap(mmu_filter, star_tool, genome_ref, end1, end2, sample, log_dir, threads, novosort)
