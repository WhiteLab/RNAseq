#!/usr/bin/python

import getpass
import subprocess
import sys

import synapseclient
from synapseclient import Project, Folder, File

from date_time import date_time

if len(sys.argv) != 4:
    sys.stderr.write('Usage:  ' + sys.argv[0] + ' {fixed fields list}{variable fields table}{container}\n')
    exit(1)
# authenticate to synapse

usr = getpass.getpass('Username: ')
password = getpass.getpass('Password: ')
sys.stderr.write(date_time() + 'Authenticating...\n')

setproxy = 'export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128; export no_proxy="rados-bionimbus-pdc.opensciencedatacloud.org"'
subprocess.call(setproxy, shell=True)
syn = synapseclient.login(usr, password)

sys.stderr.write(date_time() + 'Authentication completed.\n')
syn_id = 'syn3270015'
# create fixed lists for shared fields
fixed_fh = open(sys.argv[1], 'r')
fixed_vals = []
fixed_keys = []
for line in fixed_fh:
    line = line.rstrip('\n')
    cols = line.split('\t')
    fixed_vals.append(cols[1])
    fixed_keys.append(cols[0])
fixed_fh.close()
# open table and upload files and add fields dynamically
tbl = open(sys.argv[2], 'r')
cont = sys.argv[3]
head = next(tbl)
head = head.rstrip('\n')
lbl = head.split('\t')
src_cmd = 'unset http_proxy; unset https_proxy;. ~/.novarc;'
for line in tbl:
    line = line.rstrip('\n')
    cols = line.split('\t')
    # create new file name
    # consortium_study_group_tissueTypeAbrv_dataType_(dataSubType)*_platform_(sampleID)*
    fn = '_'.join(fixed_vals[0:3]) + '_' + '_'.join((cols[5], fixed_vals[3], cols[6], cols[0]))
    fa = fn + '_aligned.bam'
    fu = fn + '_unaligned.bam'
    try:
        swift_cmd = src_cmd + 'swift download ' + cont + ' ' + cols[-2] + ' --output ' + fa
        sys.stderr.write(date_time() + swift_cmd + '\nDownloading aligned file ' + cols[-2])
        check = subprocess.check_output(swift_cmd, shell=True, stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + 'Download for bid ' + cols[0] + ' failed. Skipping!\n')
        continue
    swift_cmd = src_cmd + 'swift download ' + cont + ' ' + cols[-1] + ' --output ' + fu
    try:
        swift_cmd = src_cmd + 'swift download ' + cont + ' ' + cols[-1] + ' --output ' + fu
        sys.stderr.write(date_time() + swift_cmd + '\nDownloading unaligned file ' + cols[-1])
        check = subprocess.check_output(swift_cmd, shell=True, stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + 'Download for bid ' + cols[0] + ' failed. Skipping!\n')
        continue
    # comment switching betwen using a local file system or ceph
    #    ln_fa='ln -s ' + cols[-2] + ' ' + fa
    #    ln_fu='ln -s ' + cols[-1] + ' ' + fu
    #    subprocess.call(ln_fa,shell=True)
    #    subprocess.call(ln_fu,shell=True)
    # sys.stderr.write('Created links ' + ln_fa + ' ' + ln_fu + '\n')
    setproxy = 'export http_proxy=http://cloud-proxy:3128; export https_proxy=http://cloud-proxy:3128; export no_proxy="rados-bionimbus-pdc.opensciencedatacloud.org"'
    subprocess.call(setproxy, shell=True)
    fao = File(fa, parent=syn_id)
    fuo = File(fu, parent=syn_id)
    for i in xrange(0, len(fixed_keys), 1):
        setattr(fao, fixed_keys[i], fixed_vals[i])
        setattr(fuo, fixed_keys[i], fixed_vals[i])
        sys.stderr.write(fixed_keys[i] + '\t' + fixed_vals[i] + '\n')
    for i in xrange(2, 8, 1):
        setattr(fao, lbl[i], cols[i])
        setattr(fuo, lbl[i], cols[i])
        sys.stderr.write(lbl[i] + '\t' + cols[i] + '\n')
    sys.stderr.write(date_time() + 'Uploading data\n')
    try:
        fao = syn.store(fao)
        fuo = syn.store(fuo)
    except:
        sys.stderr.write(date_time() + 'Upload failed for bid ' + cols[0] + '\n')
        exit(1)
    sys.stderr.write(date_time() + 'Uploading data complete for bid ' + cols[0] + '\nRemoving downloaded bam\n')
    rm_bam = 'rm ' + fa + ' ' + fu
    subprocess.call(rm_bam, shell=True)
tbl.close()
