#!/usr/bin/python
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/')
from utility.date_time import date_time
from subprocess import call
from utility.log import log


def picard_insert_size(java_tool, picard_tool, sample, log_dir):
    picard_insert_size_cmd = java_tool + " -Xmx2g -jar " + picard_tool + " CollectInsertSizeMetrics I=" + sample + ".srt.bam H=" + sample + ".insert_metrics.pdf O=" + sample + ".insert_metrics.hist  > " + log_dir + sample + ".picard.insert_size.log 2>&1"
    log(log_dir + sample + ".picard.insert_size.log", date_time() + picard_insert_size_cmd + "\n")
    call(picard_insert_size_cmd, shell=True)
    # open file and return insert size
    fh = open(sample + ".insert_metrics.hist", 'r')
    for i in xrange(0, 7, 1):
        skip = next(fh)
    stats = next(fh)
    fh.close()
    stat = stats.split('\t')

    return (stat[4], stat[5])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Picard collect insert size metrics module.  Gathers insert size metrics, run after removing BAM duplicates.')
    parser.add_argument('-j', '--java', action='store', dest='java_tool',
                        help='Java location directory, version jdk1.7.0_45 preferred')
    parser.add_argument('-p', '--picard', action='store', dest='picard_tool', help='Picard jar file location')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/project name prefix')
    parser.add_argument('-l', '--log', action='store', dest='log_dir', help='LOG directory location')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (java_tool, picard_tool, sample, log_dir) = (inputs.java_tool, inputs.picard_tool, inputs.sample, inputs.log_dir)
    picard_insert_size(java_tool, picard_tool, sample, log_dir)
