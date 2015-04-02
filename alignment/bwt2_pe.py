#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log
import subprocess

def bwt2_pe(bwt_tool,bwt_ref,end1,end2,samtools_tool,samtools_ref,sample,log_dir):
    bwt_cmd="(" + bwt_tool + " --fr -p -l 0 -X 500 -x " + bwt_ref + " -1 " + end1 + " -2 " + end2 + " | " + samtools_tool + " view -bT " + samtools_ref + " - > " + sample + ".bam) > " + log_dir + sample + ".bwt.pe.log 2>&1"
    loc=log_dir + sample + ".bwt.pe.log"
    log(loc,date_time() + bwt_cmd + "\n")
    try:
        subprocess.check_output(bwa_cmd,shell=True)
    except:
        exit(1)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Bowtie2 paired-end alignment module.  Typically run first in pipeline.')
    parser.add_argument('-b','--bwt',action='store',dest='bwt_tool', help='Location of bowtie2 alignment tool.')
    parser.add_argument('-br','--bwt_reference',action='store',dest='bwt_ref',help='Location of bwt reference file')
    parser.add_argument('-f1','--file1',action='store',dest='end1',help='First of paired-end fastq file')
    parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second of paired-end fastq file')
    parser.add_argument('-s','--samtools',action='store',dest='samtools_tool',help='Location of samtools tool.  Version 1.19 preferred.')
    parser.add_argument('-sr','--samtools_reference',action='store',dest='samtools_ref',help='Location of samtools reference')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()
    (bwt_tool,bwt_ref,end1,end2,samtools_tool,samtools_ref,sample,log_dir)=(inputs.bwt_tool,inputs.bwt_ref,inputs.end1,inputs.end2,inputs.samtools_tool,inputs.samtools_ref,inputs.sample,inputs.log_dir)
    bwt2_pe(bwt_tool,bwt_ref,end1,end2,samtools_tool,samtools_ref,sample,log_dir)
