#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log


def star(STAR, genome, end1, end2, sample, log_dir, th, sf):
    loc = log_dir + sample + "star.log"
    meta = sample.split('_')
    RGRP = "ID:" + sample + "\tLB:" + meta[0] + "\tPU:" + meta[4] + "\tSM:" + meta[0] + "\tPL:illumina"
    star_cmd = STAR + " --runMode alignReads --twopassMode Basic --outFileNamePrefix " + sample + ". --runThreadN " \
               + th + " --genomeDir " + genome + " --readFilesIn " + end1 + " " + end2 + " --readFilesCommand zcat \
               --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outFilterType BySJout \
               --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 8 \
               --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattrRGline \"" + RGRP \
               + "\""
    if sf == 'N':
        # add XS tag is input is not stranded
        star_cmd += ' --outSAMattributes NH HI AS nM XS'
    star_cmd += ' 2>> ' + loc + ' >> ' + loc + '; mv *Log* ' + log_dir

    log(loc, date_time() + star_cmd + "\n")
    check = call(star_cmd, shell=True)
    if check == 0:
        return 0
    else:
        return 1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='star paired-end alignment module.')
    parser.add_argument('-s', '--star', action='store', dest='star', help='Location of star aligner.')
    parser.add_argument('-g', '--genome', action='store', dest='genome',
                        help='Location of directory contaning genome built by STAR')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First of paired-end fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second of paired-end fastq file')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-th', '--threads', action='store', dest='th', help='Number of threads')
    parser.add_argument('-sf', '--stranded', action='store', dest='sf', help='Flag whether data is stranded')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (STAR, genome, end1, end2, sample, log_dir, th, sf) = (inputs.star, inputs.genome, inputs.end1, inputs.end2,
                                                       inputs.sample, inputs.log_dir, inputs.th, inputs.sf)
    star(STAR, genome, end1, end2, sample, log_dir, th, sf)
