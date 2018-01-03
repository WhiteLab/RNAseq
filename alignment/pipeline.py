#!/usr/bin/env python3
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')
import os
import re
from utility.date_time import date_time
import subprocess
import json
from utility.log import log
from alignment.cutadapter import cutadapter
from alignment.fastqc import fastqc
from alignment.bwt2_pe import bwt2_pe
from alignment.novosort_sort_pe import novosort_sort_pe
from alignment.picard_insert_size import picard_insert_size
from alignment.qc_bam import qc_bam
from alignment.filter_wrap import filter_wrap
from alignment.star import star
from subprocess import call
from alignment.parse_qc import parse_qc


class Pipeline:
    def __init__(self, end1, end2, json_config):
        self.json_config = json_config
        self.sf1 = end1
        self.sf2 = end2
        self.end1 = os.path.basename(self.sf1)
        self.end2 = os.path.basename(self.sf2)
        self.status = 0
        self.config_data = json.loads(open(self.json_config, 'r').read())
        s = re.match('^(\S+)_1_sequence\.txt\.gz$', self.end1)
        if s:
            self.sample = s.group(1)
        else:
            s = re.match('(^\S+)_\D*\d\.f\w*q\.gz$', self.end1)
            self.sample = s.group(1)

        HGACID = self.sample.split("_")
        self.bnid = HGACID[0]
        self.fastqc_tool = self.config_data['tools']['fastqc']
        self.java_tool = self.config_data['tools']['java']
        self.picard_tool = self.config_data['tools']['picard']
        self.novosort = self.config_data['tools']['novosort']
        self.bwt2_tool = self.config_data['tools']['bwt2']
        self.bwt2_ref = self.config_data['refs']['bwt2']
        self.samtools_tool = self.config_data['tools']['samtools']
        self.picard_tmp = 'picard_tmp'
        self.star_tool = self.config_data['tools']['star']
        self.pdxflag = self.config_data['params']['pdxflag']
        self.skip_cut = self.config_data['params']['skip_cut']
        if self.pdxflag == 'Y':
            self.mmu_filter = self.config_data['tools']['mouse_filter']
            self.mmu_star_ref = self.config_data['refs']['mmu_star']
            self.hsa_star_ref =  self.config_data['refs']['hsa_star']
        self.genome_ref = self.config_data['refs']['genome']
        self.samtools_ref =  self.config_data['refs']['samtools']
        self.gtf_ref = self.config_data['refs']['gtf']
        self.tx =  self.config_data['refs']['transcriptome']
        self.align_dir = self.config_data['refs']['align_dir']
        self.project = self.config_data['refs']['project']
        self.project_dir = self.config_data['refs']['project_dir']
        self.cwd = self.project_dir + self.project + '/' + self.align_dir + '/' + self.bnid + '/' + self.sample
        self.user = self.config_data['params']['user']
        self.group = self.config_data['params']['group']
        self.threads = self.config_data['params']['threads']
        self.ram = self.config_data['params']['ram']
        self.sf = self.config_data['params']['stranded']
        self.bam_dir = 'BAMS/'
        self.qc_dir = 'QC/'
        self.log_dir = 'LOGS/'
        self.star_dir = 'STAR_OUT/'
        self.fq_trimmed = 'TRIMMED_FQ/'
        self.loc = self.log_dir + self.sample + '.pipe.log'
        self.pipeline()

    def pipeline(self):
        # temp line to source environment variables until compute is restarted
        src_env = '. /etc/environment'
        call(src_env, shell=True)

        # create working directory
        if not os.path.isdir(self.cwd):
            mk_cwd = 'mkdir -p ' + self.cwd
            sys.stderr.write(date_time() + 'Creating working directory ' + mk_cwd + '\n')
            call(mk_cwd, shell=True)
        os.chdir(self.cwd)
        if not os.path.isdir(self.log_dir):
            mk_log_dir = 'mkdir ' + self.log_dir
            call(mk_log_dir, shell=True)
            log(self.loc, date_time() + 'Made log directory ' + self.log_dir + "\n")
        if not os.path.isdir(self.fq_trimmed):
            mk_fq_dir = 'mkdir ' + self.fq_trimmed
            call(mk_fq_dir, shell=True)
            log(self.loc, date_time() + 'Made fastq trimmed directory ' + self.fq_trimmed + "\n")
        os.chdir(self.fq_trimmed)
        mv_fq = 'mv ../' + self.log_dir + ' .'
        call(mv_fq, shell=True)
        log(self.loc, date_time() + 'Changed into ' + self.fq_trimmed + " and moved log directory there\n")
        if not os.path.isdir(self.star_dir):
            mk_star_dir = 'mkdir ' + self.star_dir
            call(mk_star_dir, shell=True)
            log(self.loc, date_time() + 'Made star output directory ' + self.star_dir + "\n")
        if not os.path.isdir(self.bam_dir):
            mk_bam_dir = 'mkdir ' + self.bam_dir
            call(mk_bam_dir, shell=True)
            log(self.loc, date_time() + 'Made bam output directory ' + self.bam_dir + "\n")
        if not os.path.isdir(self.qc_dir):
            mk_qc_dir = 'mkdir ' + self.qc_dir
            call(mk_qc_dir, shell=True)
            log(self.loc, date_time() + 'Made qc directory ' + self.qc_dir + "\n")
        log(self.loc,
            date_time() + "Starting alignment qc for paired end sample files " + self.end1 + " and " + self.end2 + "\n")

        # remove adapters
        if self.skip_cut == 'N':
            check = cutadapter(self.sample, self.sf1, self.sf2, self.json_config)
            if check != 0:
                log(self.loc, date_time() + 'cutadapt failure for ' + self.sample + '\n')
                exit(1)

            # start fastqc, will run while insert size being calculated

            end_ss1 = self.sample + '_1.subset.fastq'
            end_ss2 = self.sample + '_2.subset.fastq'
            subset = self.sample + '_subset'

            ss_cmd = 'gunzip -c ' + self.end1 + ' | head -n 4000000 > ' + end_ss1
            subprocess.call(ss_cmd, shell=True)
            ss_cmd = 'gunzip -c ' + self.end2 + ' | head -n 4000000 > ' + end_ss2
            subprocess.call(ss_cmd, shell=True)
            # check certain key processes

            check = bwt2_pe(self.bwt2_tool, self.tx, end_ss1, end_ss2, self.samtools_tool, self.samtools_ref, subset,
                            self.threads, self.log_dir)
            if check != 0:
                log(self.loc, date_time() + 'Bowtie2 failure for ' + self.sample + '\n')
                self.status = 1
                exit(1)
            check = novosort_sort_pe(self.novosort, subset, self.log_dir, self.threads,
                                     self.ram, 'coord')  # rest won't run until completed
            if check != 0:
                log(self.loc, date_time() + 'novosort sort failure for ' + self.sample + '\n')
                self.status = 1
                exit(1)
            (x, s) = picard_insert_size(self.java_tool, self.picard_tool, subset, self.log_dir)
            log(self.loc, date_time() + 'Running qc on fastq file\n')
            fastqc(self.fastqc_tool, self.sample, self.end1, self.end2, self.threads)
        if self.pdxflag == 'Y':
            log(self.loc, date_time() + 'Aligning and filtering reads for mouse contamination')
            check = filter_wrap(self.mmu_filter, self.star_tool, self.mmu_star_ref, self.end1, self.end2,
                                self.sample, self.log_dir, self.threads, self.novosort)
            if check != 0:
                log(self.loc, date_time() + 'Read filter failure for ' + self.sample + '\n')
                exit(1)
            log(self.loc, date_time() + 'Performing star alignment ' + self.sample + '\n')
            check = star(self.star_tool, self.hsa_star_ref, self.end1, self.end2, self.sample, self.log_dir,
                         self.threads, self.sf)
        else:
            log(self.loc, date_time() + 'Starting star align\n')
            check = star(self.star_tool, self.genome_ref, self.end1, self.end2, self.sample, self.log_dir, self.threads,
                     self.sf)

        if check != 0:
            log(self.loc, date_time() + 'star alignment failure for ' + self.sample + '\n')
            self.status = 1
            exit(1)
        # run QC on bams
        check = qc_bam(self.sample, self.json_config)
        if check != 0:
            log(self.loc, date_time() + 'bam qc process failure for ' + self.sample + '\n')
            self.status = 1
            exit(1)

        check = parse_qc(self.json_config, self.sample)
        if check != 0:
            log(self.loc, date_time() + 'qc summary failure for ' + self.sample + '\n')
            self.status = 1
            exit(1)
        # move outputs to correct directories
        log(self.loc, date_time() + 'Organizing outputs\n')
        mv_bams = 'mv *Aligned*.bam  ' + self.bam_dir
        call(mv_bams, shell=True)
        mv_star = 'mv  *.tab ' + self.star_dir
        call(mv_star, shell=True)
        mv_sub = 'mv *subset.insert* *.txt *.pdf *.json ' + self.qc_dir + '; cp ' + self.json_config + ' ' + self.qc_dir
        call(mv_sub, shell=True)
        # mv subdirectories to right place
        mv_dir = 'mv ' + ' '.join((self.bam_dir, self.log_dir, self.qc_dir, self.star_dir)) + ' ../'
        call(mv_dir, shell=True)
        rm_tmp = 'rm -rf *STAR* *subset*'
        call(rm_tmp, shell=True)
        set_acl = 'chown -R ' + self.user + ':' + self.group + ' ./;'
        sys.stderr.write(date_time() + 'Setting acls for current directory ' + set_acl + '\n')
        call(set_acl, shell=True)
        os.chdir('../../')
        sys.stderr.write(date_time() + 'Moving subdirs out of ' + self.sample + '\n')
        mv_dirs = 'mv ' + self.sample + '/* .'
        call(mv_dirs, shell=True)
        rm_lane = 'rmdir ' + self.sample
        call(rm_lane, shell=True)

        sys.stderr.write(date_time() + 'Pipeline complete for ' + self.sample + '\n')
        self.status = 0


def main():
    import argparse
    parser = argparse.ArgumentParser(description='RNA alignment paired-end QC pipeline')
    parser.add_argument('-f1', '--file1', action='store', dest='end1', help='First fastq file')
    parser.add_argument('-f2', '--file2', action='store', dest='end2', help='Second fastq file')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file containing tool and reference locations')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()

    end1 = inputs.end1
    end2 = inputs.end2
    config_file = inputs.config_file
    Pipeline(end1, end2, config_file)


if __name__ == "__main__":
    main()
