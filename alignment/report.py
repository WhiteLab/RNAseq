#!/usr/bin/python
import sys

sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
import os
import re
from utility.date_time import date_time
from utility.log import log


def report(sample, ref_gtf, tx_gtf):
    # casual logging - look for a LOGS directory, otherwise assume current dir
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.report.log'
    log(loc, date_time() + "Loading reference transcripts from " + ref_gtf + "\n")
    ref = {}
    rh = open(ref_gtf, 'r')
    for line in rh:
        if line[0:2] == '##':
            continue
        attr = line.split('\t')
        (bio, tx_id, name) = ('', '', '')
        # might have mixed type between gencode and ensenmbl, test
        #        pdb.set_trace()
        if re.search('gene_biotype "(\S+)";', attr[-1]):
            m = re.search('gene_biotype "(\S+)";.*transcript_id "(\S+)"; transcript_name "(\S+)";', attr[-1])

            try:
                (bio, tx_id, name) = (m.group(1), m.group(2), m.group(3))
            except:
                sys.stderr.write('Pattern failed for line ' + line)
                m = re.search('gene_biotype "(\S+)";.*transcript_id "(\S+)"; transcript_name (\S+);', attr[-1])
            #                exit(1)
        else:
            m = re.search('transcript_id "(\S+)"; gene_type "(\S+)";.* gene_name "(\S+)";', attr[-1])
            (bio, tx_id, name) = (m.group(2), m.group(1), m.group(3))

        if tx_id not in ref:
            ref[tx_id] = {}
            ref[tx_id]['type'] = bio
            ref[tx_id]['name'] = name
    rh.close()
    log(loc,
        date_time() + 'Completed loading annotation, parsing and annotating cufflinks output file ' + tx_gtf + '\n')
    rpt = open(sample + '.report.txt', 'w')
    rpt.write(
        'chr\tstrand\tstart\tend\tscore\tgeneid\ttxid\tgene_sym\tgene_biotype\tfpkm\tconflo\tconfhi\tcoverage\tread_support\n')
    th = open(tx_gtf, 'r')
    flist = ('FPKM', 'conf_lo', 'conf_hi', 'cov', 'full_read_support')
    chrom = ''
    for line in th:
        line = line.rstrip('\n')
        anno = line.split('\t')
        if anno[0] != chrom:
            log(loc, date_time() + 'Processing results for chromosome ' + chrom + '\n')
            chrom = anno[0]
        if anno[2] == 'transcript':
            m = re.findall('(\w+) \"(\S+)\";', anno[-1])
            fields = {}
            for i in xrange(0, len(m), 1):
                fields[m[i][0]] = m[i][1]
            # only print values where fpkm > 0
            if float(fields['FPKM']) > 0.0:
                rpt.write('\t'.join([anno[0], anno[6], anno[3], anno[4], anno[5]]) + '\t')
                #                pdb.set_trace()
                rpt.write(fields['gene_id'] + '\t' + fields['transcript_id'] + '\t')
                if fields['transcript_id'] in ref:
                    rpt.write(ref[fields['transcript_id']]['name'] + '\t' + ref[fields['transcript_id']]['type'])
                else:
                    rpt.write('NA\tNA')
                for i in xrange(0, len(flist), 1):
                    rpt.write('\t' + fields[flist[i]])
                rpt.write('\n')
    th.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Transcript annotation summart report.  Parses cufflinks trnascripts.gtf and adds gene symbols and transcript biotype')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample/location name prefix')
    parser.add_argument('-r', '--reference', action='store', dest='ref_gtf',
                        help='Reference gtf file used to annotate - preferably GENCODE with gene name and gene type')
    parser.add_argument('-c', '--cufflinks', action='store', dest='tx_gtf', help='Cufflinks transcript output file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (sample, ref_gtf, tx_gtf) = (inputs.sample, inputs.ref_gtf, inputs.tx_gtf)
    report(sample, ref_gtf, tx_gtf)
