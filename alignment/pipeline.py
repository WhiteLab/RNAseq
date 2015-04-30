#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/alignment')
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import os
import re
from date_time import date_time
from fastqc import fastqc
from bwt2_pe import bwt2_pe
#from picard_sort_pe import picard_sort_pe
from novosort_sort_pe import novosort_sort_pe
from picard_insert_size import picard_insert_size
from tophat import tophat
from cufflinks import cufflinks
from subprocess import call
import subprocess
import json
from log import log
import pdb

class Pipeline():
    def __init__(self,end1,end2,seqtype,json_config,ref_mnt):
        self.json_config = json_config
        self.end1=end1
        self.end2=end2
        self.seqtype=seqtype
        self.status=0
        self.ref_mnt=ref_mnt
        self.parse_config()

    def parse_config(self):
        self.config_data = json.loads(open(self.json_config, 'r').read())
        s=re.match('^(\S+)_1_sequence\.txt\.gz$',self.end1)
        if s:
            self.sample=s.group(1)
        else:
            s=re.match('(^\S+)_\D*\d\.f\w*q\.gz$',self.end1)
            self.sample=s.group(1)
        self.loc='LOGS/' + self.sample + '.pipe.log'
        HGACID=self.sample.split("_")
        self.bid=HGACID[0]
        self.fastx_tool=self.config_data['tools']['fastx']
        self.bwt2_tool=self.config_data['tools']['bwt2']
        self.bwa_ref=self.ref_mnt + '/' + self.config_data['refs']['bwa']
        self.samtools_tool=self.config_data['tools']['samtools']
        self.samtools_ref=self.ref_mnt + '/' + self.config_data['refs']['samtools']
        self.java_tool=self.config_data['tools']['java']
        self.picard_tool=self.config_data['tools']['picard']
        self.novosort=self.config_data['tools']['novosort']
        self.picard_tmp='picard_tmp'
        self.tophat=self.config_data['tools']['tophat']
        self.cufflinks=self.config_data['tools']['cufflinks']
        self.htseq_count=self.config_data['tools']['htseq-count']
        self.bedtools2_tool=self.config_data['tools']['bedtools']
        self.bed_ref=self.ref_mnt + '/' + self.config_data['refs'][self.seqtype]
        self.obj=self.config_data['refs']['obj']
        self.cont=self.config_data['refs']['cont']
        self.pipeline()

    def pipeline(self):
        log_dir='LOGS/'
        if os.path.isdir(log_dir) == False:
            mk_log_dir='mkdir ' + log_dir
            call(mk_log_dir,shell=True)
            log(self.loc,date_time() + 'Made log directory ' + log_dir + "\n")
        # create BAM and QC directories if they don't exist already
        bam_dir='BAM/'
        qc_dir='QC/'
        if os.path.isdir(bam_dir) == False:
            mk_bam_dir='mkdir ' + bam_dir
            call(mk_bam_dir,shell=True)
            log(self.loc,date_time() + 'Made bam directory ' + bam_dir + "\n")
        if os.path.isdir(qc_dir) == False:
            mk_qc_dir='mkdir ' + qc_dir
            call(mk_qc_dir,shell=True)
            log(self.loc,date_time() + 'Made qc directory ' + qc_dir + "\n")
        log(self.loc,date_time() + "Starting alignment qc for paired end sample files " + self.end1 + " and " + self.end2 + "\n")
        #inputs
        
        SAMPLES={}
        SAMPLES[self.sample]={}
        SAMPLES[self.sample]['f1']=self.end1
        SAMPLES[self.sample]['f2']=self.end2
        RGRP="@RG\\tID:" + self.sample + "\\tLB:" + self.bid + "\\tSM:" + self.bid + "\\tPL:illumina"
        
        # use subset of fastq files to get insert size estimate
        end_ss1=sample + '_1.subset.fastq'
        end_ss2=sample + '_2.subset.fastq'
        ss_cmd='gunzip -c ' + end1 + ' > ' + end_ss1
        subprocess.call(ss_cmd,shell=True)
        ss_cmd='gunzip -c ' + end2 + ' > ' + end_ss2
        subprocess.call(ss_cmd,shell=True)
        # check certain key processes
        check=bwt2_pe(self.bwa_tool,RGRP,self.bwa_ref,end_ss1,end_ss2,self.samtools_tool,self.samtools_ref,self.sample,log_dir) # rest won't run until completed
        if(check != 0):
            log(self.loc,date_time() + 'Bowtie2 failure for ' + self.sample + '\n')
            exit(1)
        check=novosort_sort_pe(self.novosort,self.sample,log_dir) # rest won't run until completed
        if(check != 0):
            log(self.loc,date_time() + 'novosort sort failure for ' + self.sample + '\n')
            exit(1)
        # start fastqc, will run while insert size being calculated
        log(self.loc,date_time() + 'Running qc on fastq file\n')
        fastqc(self.fastqc_tool,self.sample,self.end1,self.end2) # flag determines whether to run independently or hold up the rest of the pipe until completion
        log(self.loc,date_time() + 'Calculating insert sizes\n')
        (x,s)=picard_insert_size(self.java_tool,self.picard_tool,self.sample,log_dir) # get insert size metrics, use for tophat. 
        log(self.loc,date_time() + 'Performing tophat alignment ' + self.sample + '\n')
        tophat(self.tophat_tool,self.ens_ref,self.bwt2_ref,self.end1,self,end2,x,s,self.sample,log_dir)

def main():
    import argparse
    parser=argparse.ArgumentParser(description='RNA alignment paired-end QC pipeline')
    parser.add_argument('-f1','--file1',action='store',dest='end1',help='First fastq file')
    parser.add_argument('-f2','--file2',action='store',dest='end2',help='Second fastq file')
    parser.add_argument('-t','--seqtype',action='store',dest='seqtype',help='Type of sequencing peformed.  Likely choices are genome, exome, and capture')
    parser.add_argument('-j','--json',action='store',dest='config_file',help='JSON config file containing tool and reference locations')
    parser.add_argument('-m','--mount',action='store',dest='ref_mnt',help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    inputs=parser.parse_args()

    end1=inputs.end1
    end2=inputs.end2
    seqtype=inputs.seqtype
    config_file=inputs.config_file
    ref_mnt=inputs.ref_mnt
    Pipeline(end1,end2,seqtype,config_file,ref_mnt)
if __name__ == "__main__":
    main()
