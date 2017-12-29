#!/usr/bin/env python3
import sys
sys.path.append('/cephfs/users/mbrown/RNAseq')
from utility.date_time import date_time
import subprocess
from utility.log import log
import os


def novosort_sort_pe(novosort, sample, log_dir, t, mem, temp, stype):
    samp_root = os.path.basename(sample)
    novosort_sort_pe_cmd = 'mkdir ' + temp + ';' + novosort + " --threads " + t + " --ram " + mem \
                           + "G --tmpdir  " + temp + " --output " + sample + ".srt.bam --index  " + sample + ".bam > " \
                           + log_dir + samp_root + ".novosort.sort.pe.log 2>&1"
    if stype == 'name':
        novosort_sort_pe_cmd = 'mkdir ' + temp + ';' + novosort + " --threads " + t + " --ram " + mem \
                               + "G --tmpdir  " + temp + " --output " + sample + ".nsrt.bam -n  " + sample + ".bam > " \
                               + log_dir + samp_root + ".novosort.sort.pe.log 2>&1"
    log(log_dir + samp_root + ".novosort.sort.pe.log", date_time() + novosort_sort_pe_cmd + "\n")
    f = 0
    try:
        f = subprocess.call(novosort_sort_pe_cmd, shell=True)
        #rm_tmp = 'rm -rf novosort_tmp'
        #subprocess.call(rm_tmp, shell=True)
    except:
        log(log_dir + sample + ".novosort.sort.pe.log", 'novosort sort failed for sample ' + sample + '\n')
        exit(1)
    return f


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='novosort tool to sort BAM module.')
    parser.add_argument('-n', '--novosort', action='store', dest='novosort', help='novosort binary location')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')
    parser.add_argument('-t', '--threads', action='store', dest='t', help='Number of threads')
    parser.add_argument('-m', '--memory', action='store', dest='mem', help='Memory - in GB')
    parser.add_argument('-d', '--temp_dir', action='store', dest='temp', help='temp dir for sort')
    parser.add_argument('-y', '--sort_type', action='store', dest='stype', help='type of sort, name or coord')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (novosort, sample, log_dir, t, mem, temp, stype) = (inputs.novosort, inputs.sample, inputs.log_dir, inputs.t, inputs.mem,
                                                 inputs.temp, inputs.stype)
    novosort_sort_pe(novosort, sample, log_dir, t, mem, temp, stype)
