#!/usr/bin/python
# written by Miguel Brown 2015-Feb-23. Wrapper script to loop through sequencing files and use pipeline

import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
import os
import re
import argparse
from job_manager import job_manager
from date_time import date_time
import subprocess
import pdb

def batch_qc(fn,cont,obj,t):
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    inputs=parser.parse_args()
    fh=open(inputs.fn,'r')
    src_cmd='. ~/.novarc;'
    jobs=[]
    for line in fh:
        line=line.rstrip('\n')    
        # All files for current bid to be stored in cwd
        swift_cmd=src_cmd + 'swift list ' + cont + ' --prefix ' + obj + '/' + line
        sys.stderr.write(date_time() + 'Checking for sequence files for sample ' + line + '\n' + swift_cmd + '\n')
        try:
            contents=subprocess.check_output(swift_cmd,shell=True)
            if len(contents) < len(line):
                sys.stderr.write(date_time() + 'Can\'t find sequencing files for ' + line + ' skipping!\n')
                continue
        except:
            sys.stderr.write(date_time() + 'Can\'t find sequencing files for ' + line + ' skipping!\n')
            continue
        seqfile=re.findall('(\S+[sequence|f*q]*\.gz)',contents)
        sf1=seqfile[0]
        end1=os.path.basename(sf1)
        sf2=seqfile[1]
        end2=os.path.basename(sf2)
        swift_cmd=src_cmd + "swift download " + cont + " --skip-identical --prefix " + obj + '/' + line
        link_cmd='ln -s ' + sf1 + ' .;ln -s ' + sf2
        fastqc_cmd='mkdir -p PREQC/' + line + '; fastqc -t 2 -o PREQC/' + line + ' ' + sf1 + ' ' + sf2
        upload_cmd=src_cmd + 'swift upload ' + cont + ' PREQC/' + line
        cleanup_cmd='rm -rf RAW/' + line + ' PREQC/' + line + ' ' + end1 + ' ' + end2
        jobs.append(';'.join([swift_cmd,link_cmd,fastqc_cmd,upload_cmd,cleanup_cmd]))
    sys.stderr.write(date_time() + 'Job list created, running jobs!\n')
    job_manager(jobs,t)
    return 0

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Pipeline wrapper script to process multiple paired end set serially.')
    parser.add_argument('-f','--file',action='store',dest='fn',help='File with bionimbus ID list')
    parser.add_argument('-c','--container',action='store',dest='cont',help='Name of swift container')
    parser.add_argument('-o','--object',action='store',dest='obj',help='Swift object prefix, like RAW/')
    parser.add_argument('-t','--threads',action='store',dest='t',help='Number of simultaneous threads to keep going')
    inputs=parser.parse_args()
    (fn,cont,obj,t)=(inputs.fn,inputs.cont,inputs.obj,inputs.t)
    batch_qc(fn,cont,obj,t)
