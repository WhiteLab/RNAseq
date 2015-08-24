#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log
import subprocess

def star(STAR,genome,end1,end2,sample,log_dir,th):
    loc=log_dir + sample + "star.log"
    meta=sample.split('_')
    epoch=150409
    star_cmd=STAR + " --runMode alignReads --twopassMode --outFileNamePrefix " + sample + " --runThreadN " + th + " --genomeDir " + genome + " --readFilesIn " + end1 + " " + end2 + " --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 8 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 2>> " + loc
    
    log(loc,date_time() + star_cmd + "\n")
    call(star_cmd,shell=True)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='star paired-end alignment module.')
    parser.add_argument('-s','--star',action='store',dest='star', help='Location of star aligner.')
    parser.add_argument('-g','--genome',action='store',dest='genome',help='Location of directory contaning genome built by STAR')
    parser.add_argument('-f1','--file1',action='store',dest='end1',help='First of paired-end fastq file')
    parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second of paired-end fastq file')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')
    parser.add_argument('-th','--threads',action='store',dest='th',help='Number of threads')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    inputs=parser.parse_args()
    (STAR,sam,genome,end1,end2,sample,log_dir,th)=(inputs.star,inputs.genome,inputs.end1,inputs.end2,inputs.sample,inputs.log_dir,inputs.th)
    star(STAR,genome,end1,end2,sample,log_dir,th)
