#!/usr/bin/python

import subprocess
import sys
from subprocess import check_output

from date_time import date_time


def download_from_swift(cont, obj):
    obj = obj.rstrip('\n')
    src_cmd = ". /home/ubuntu/.novarc;"
    swift_cmd = src_cmd + "swift download " + cont + " --skip-identical --prefix " + obj
    sys.stderr.write(date_time() + swift_cmd + "\n")
    try:
        check = check_output(swift_cmd, shell=True, stderr=subprocess.PIPE)
    except:
        sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
        exit(1)
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Simple download module to get files from swift.  Can use prefix or whole object name')
    parser.add_argument('-c', '--container', action='store', dest='project', help='Swift container, i.e. PANCAN')
    parser.add_argument('-o', '--object', action='store', dest='align_dir',
                        help='Swift object name/prefix, i.e. RAW/2015-1234')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (cont, obj) = (inputs.cont, inputs.obj)
    download_from_swift(cont, obj)
