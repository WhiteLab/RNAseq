#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from subprocess import call
from log import log
import subprocess

def cufflinks(cufflinks_tool,ens_ref,sample,log_dir):
    loc=log_dir + sample + ".cufflinks.log"
    cufflinks_cmd=cufflinks_tool + " --GTF " + ens_ref + " -p 8 --library-type fr-firststrand --frag-bias-correct --multi-read-correct --upper-quartile-norm --pre-mrna-fraction -o " + sample + " 2>> " + loc
    log(loc,date_time() + cufflinks_cmd + "\n")
    try:
        subprocess.check_output(cufflinks_cmd,shell=True)
    except:
        exit(1)
    return 0

if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description='tophat paired-end alignment module.  Typically run first in pipeline.')
    parser.add_argument('-c','--cufflinks',action='store',dest='cufflinks_tool', help='Location of cufflinks tool.')
    parser.add_argument('-e','--ensembl_reference',action='store',dest='ens_ref',help='Location of ensembl reference file')
    parser.add_argument('-sa','--sample',action='store',dest='sample',help='Sample/project name prefix')
    parser.add_argument('-l','--log',action='store',dest='log_dir',help='LOG directory location')

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    inputs=parser.parse_args()
    (cufflinks_tool,ens_ref,sample,log_dir)=(inputs.cufflinks_tool,inputs.ens_ref,inputs.sample,inputs.log_dir)
    cufflinks_tool(cufflinks_tool,ens_ref,sample,log_dir)
