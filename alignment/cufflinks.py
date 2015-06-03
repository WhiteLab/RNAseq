#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log
import subprocess

def cufflinks(cufflinks_tool,ens_ref,genome,sample,log_dir,t):
    loc=log_dir + sample + ".cufflinks.log"
    cufflinks_cmd=cufflinks_tool + " " + sample + "/accepted_hits.bam -g " + ens_ref + " -p " + t + " --library-type fr-secondstrand -b " + genome + " -u --upper-quartile-norm --pre-mrna-fraction -o " + sample + " 2>> " + loc
    log(loc,date_time() + cufflinks_cmd + "\n")
    try:
        subprocess.check_output(cufflinks_cmd,shell=True)
    except:
        exit(1)
    
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='Annotation and fpkm estimation.  Run after tophat.')
    parser.add_argument('-c','--cufflinks',action='store',dest='cufflinks_tool', help='Location of cufflinks tool.')
    parser.add_argument('-e','--ensembl_reference',action='store',dest='ens_ref',help='Location of ensembl reference file')
    parser.add_argument('-g','--genome',action='store',dest='genome',help='Location of genome reference file')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')
    parser.add_argument('-t','--threads',action='store',dest='t',help='Number of threads')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    inputs=parser.parse_args()
    (cufflinks_tool,ens_ref,genome,sample,log_dir,t)=(inputs.cufflinks_tool,inputs.ens_ref,inputs.genome,inputs.sample,inputs.log_dir,inputs.t)
    cufflinks(cufflinks_tool,ens_ref,genome,sample,log_dir,t)
