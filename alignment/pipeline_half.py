#!/usr/bin/python
import sys

import os
import re
from utility.date_time import date_time
import json
from utility.log import log
from cutadapter import cutadapter
from fastqc_mod import fastqc
from utility.upload_to_swift import upload_to_swift
from subprocess import call


class Pipeline():
    def __init__(self, end1, end2, json_config, ref_mnt):
        self.json_config = json_config
        self.end1 = end1
        self.end2 = end2
        self.status = 0
        self.ref_mnt = ref_mnt
        self.parse_config()

    def parse_config(self):
        self.config_data = json.loads(open(self.json_config, 'r').read())
        s = re.match('^(\S+)_1_sequence\.txt\.gz$', self.end1)
        if s:
            self.sample = s.group(1)
        else:
            s = re.match('(^\S+)_\D*\d\.f\w*q\.gz$', self.end1)
            self.sample = s.group(1)
        self.loc = 'LOGS/' + self.sample + '.pipe.log'
        HGACID = self.sample.split("_")
        self.bid = HGACID[0]
        self.fastqc_tool = self.config_data['tools']['fastqc']
        self.java_tool = self.config_data['tools']['java']
        self.obj = self.config_data['refs']['obj']
        self.cont = self.config_data['refs']['cont']
        self.threads = self.config_data['params']['threads']
        self.ram = self.config_data['params']['ram']
        self.pipeline()

    def pipeline(self):
        log_dir = 'LOGS/'
        if os.path.isdir(log_dir) == False:
            mk_log_dir = 'mkdir ' + log_dir
            call(mk_log_dir, shell=True)
            log(self.loc, date_time() + 'Made log directory ' + log_dir + "\n")
        fq_dir = 'TRIMMED_FQ'
        if os.path.isdir(fq_dir) == False:
            mk_fq_dir = 'mkdir ' + fq_dir
            call(mk_fq_dir, shell=True)
            log(self.loc, date_time() + 'Made fastq trimmed directory ' + fq_dir + "\n")
        os.chdir(fq_dir)
        mv_fq = 'mv ../LOGS .'
        call(mv_fq, shell=True)
        log(self.loc, date_time() + 'Changed into ' + fq_dir + " and moved log directory there\n")
        # create QC directories if they don't exist already
        qc_dir = 'QC/'
        if os.path.isdir(qc_dir) == False:
            mk_qc_dir = 'mkdir ' + qc_dir
            call(mk_qc_dir, shell=True)
            log(self.loc, date_time() + 'Made qc directory ' + qc_dir + "\n")
        log(self.loc,
            date_time() + "Starting alignment qc for paired end sample files " + self.end1 + " and " + self.end2 + "\n")
        # inputs

        SAMPLES = {}
        SAMPLES[self.sample] = {}
        SAMPLES[self.sample]['f1'] = self.end1
        SAMPLES[self.sample]['f2'] = self.end2
        # remove adapters

        check = cutadapter(self.sample, self.end1, self.end2, self.json_config)
        if (check != 0):
            log(self.loc, date_time() + 'cutadapt failure for ' + self.sample + '\n')
            exit(1)

        # start fastqc, will run while insert size being calculated

        log(self.loc, date_time() + 'Running qc on fastq file\n')
        fastqc(self.fastqc_tool, self.sample, self.end1, self.end2, self.threads)
        # move outputs to correct directories and upload
        log(self.loc, date_time() + 'Organizing outputs\n')
        rm_fq = 'rm ../' + self.end1 + ' ../' + self.end2
        call(rm_fq, shell=True)
        log(self.loc, date_time() + rm_fq + '\n')
        mv_results = 'mv ' + qc_dir + ' ' + log_dir + ' ../'
        call(mv_results, shell=True)
        os.chdir('../../')
        sys.stderr.write(date_time() + 'Uploading results for ' + self.sample + '\n')
        check = upload_to_swift(self.cont, self.obj)
        if (check != 0):
            sys, stderr, write(date_time() + 'Upload failure for ' + self.sample + '\n')
            self.status = 1
            exit(1)
        sys.stderr.write(date_time() + 'Pipeline complete for ' + self.sample + '\n')
        self.status = 0


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Just a QC pipeline running cutadapt and fastqc')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second fastq file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')
    parser.add_argument('-m', '--mount', action='store', dest='ref_mnt',
                        help='Drive mount location.  Example would be /mnt/cinder/REFS_XXX')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    end1 = inputs.end1
    end2 = inputs.end2
    config_file = inputs.config_file
    ref_mnt = inputs.ref_mnt
    Pipeline(end1, end2, config_file, ref_mnt)


if __name__ == "__main__":
    main()
