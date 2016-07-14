#!/usr/bin/env python

import json
import os
import re
import sys
import time
import pdb
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
from log import log


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['params']['r1adapt'], config_data['params']['r2adapt'], config_data['params']['minlen'], \
           config_data['params']['mqual'], config_data['params']['strand']


def skip_lines(fh, stop):
    for i in xrange(0, stop, 1):
        next(fh)


def parseFASTQC(FASTQC, loc):
    try:
        fh = open(FASTQC, 'r')
        skip_lines(fh, 8)
        len_range = next(fh)
        info = len_range.rstrip('\n').split('-')
        return info[1].rstrip('\t')
    except:
        log(loc, date_time() + 'Unable to open/process file ' + FASTQC)
        exit(1)


# dumb helper function to remove parans and retun only last part of line
def process_parens(cur):
    info = cur.rstrip('\n').split()
    data = info[-1].replace('(', '')
    data = data.replace(')', '')
    return data


def parseCUTADAPT(CUTADAPT, loc):
    try:
        pdb.set_trace()
        fh = open(CUTADAPT, 'r')
        flag = 0
        stats = []
        while flag == 0:
            cur = next(fh)
            if re.search('Total read', cur):
                # total read pairs
                stats.append(process_parens(cur))
                cur = next(fh)
                # r1a pct
                stats.append(process_parens(cur))
                cur = next(fh)
                # r2a pct
                stats.append(process_parens(cur))
                cur = next(fh)
                # too short
                stats.append(process_parens(cur))
                cur = next(fh)
                # too rp pass
                stats.append(process_parens(cur))
                next(fh)
                flag = 1
        tot_bp_line = next(fh)
        info = tot_bp_line.split()
        tot_bp = int(info[-2].replace(',', ''))
        # total bp
        stats.append(str(tot_bp))
        next(fh)
        next(fh)
        next(fh)
        # calculate trimmed base pers per read as a pct
        r1_qt_line = next(fh)
        info = r1_qt_line.split()
        r1_pct = round(float(info[-2].replace(',', ''))/tot_bp * 100)
        #r1 trimmed
        stats.append(str(r1_pct) + '%')

        r2_qt_line = next(fh)
        info = r2_qt_line.split()
        r2_pct = round(float(info[-2].replace(',', ''))/tot_bp * 100)
        # r2 trimmed
        stats.append(str(r2_pct) + '%')
        # total written
        tw = next(fh)
        stats.append(process_parens(tw))

        return stats
    except:
        log(loc, date_time() + 'Unable to open/process file ' + CUTADAPT + '\n')
        exit(1)
    #return tot_pairs, r1a_pct, r2a_pct, short, rp_pass, tot_bp, r1_trim, r2_trim, bp_pass


def parseINS(INS, loc):
    try:
        fh = open(INS, 'r')
        skip_lines(fh, 8)
        line = next(fh)
        line = line.rstrip('\n')
        stats = line.split('\t')
        fh.close()
        return stats[0], stats[1], stats[4], stats[5]
    except:
        log(loc, date_time() + 'Unable to open/process file ' + INS + '\n')
        exit(1)


def processSTAR(line):
    info = line.rstrip('\n').split('\t')
    return info[-1]


def parseSTAR(STAR, loc):
    try:
        fh = open(STAR, 'r')
        stats = []
        skip_lines(fh, 6)
        num_rds = next(fh)
        num_rds = processSTAR(num_rds)
        stats.append(num_rds)
        skip_lines(fh, 4)
        uniq = next(fh)
        uniq = processSTAR(uniq)
        stats.append(uniq)
        next(fh)
        sjt = next(fh)
        sjt = processSTAR(sjt)
        stats.append(sjt)
        skip_lines(fh, 4)
        nsj = next(fh)
        nsj = processSTAR(nsj)
        stats.append(nsj)
        mm = next(fh)
        mm = processSTAR(mm)
        stats.append(mm)
        delrate = next(fh)
        delrate = processSTAR(delrate)
        stats.append(delrate)
        next(fh)
        ins = next(fh)
        ins = processSTAR(ins)
        stats.append(ins)
        skip_lines(fh, 4)
        mml = next(fh)
        mml = processSTAR(mml)
        stats.append(mml)
        next(fh)
        mmml = next(fh)
        mmml = processSTAR(mmml)
        stats.append(mmml)
        next(fh)
        unmap = next(fh)
        unmap = processSTAR(unmap)
        unmap2 = next(fh)
        unmap2 = processSTAR(unmap2)
        unmap3 = next(fh)
        unmap3 = processSTAR(unmap3)
        unmap_tot = (float(unmap.rstrip('%')) + float(unmap2.rstrip('%')) + float(unmap3.rstrip('%')))
        stats.append(str(unmap_tot) + '%')
        return stats
    except:
        log(loc, date_time() + 'Unable to open/process file ' + STAR + '\n')
        exit(1)


def parsePICARD(PICARD, loc):
    try:
        fh = open(PICARD, 'r')
        skip_lines(fh, 8)
        keys = next(fh)
        keys = keys.rstrip('\n').split('\t')
        vals = next(fh)
        vals = vals.rstrip('\n').split('\t')
        qc_dict = {}
        for i in xrange(0, len(keys), 1):
            qc_dict[keys[i]] = vals[i]
        return qc_dict
    except:
        log(loc, date_time() + 'Unable to open/process file ' + PICARD + '\n')
        exit(1)


def parse_qc(config_file, sample):
    RG = sample.split('_')
    (r1, r2, minlen, minqual, strand) = parse_config(config_file=config_file)
    log_dir = './'
    if os.path.isdir('LOGS'):
        log_dir = 'LOGS/'
    loc = log_dir + sample + '.qc_stats.log'
    # set up input files and variables - .hist, .flagstats, .metrics, .qs
    cutadapt = log_dir + sample + '.cutadapt.log'
    star = sample + '.Log.final.out'
    insert = sample + '_subset.insert_metrics.hist'
    fastqc = 'QC/' + sample + '_1_sequence_fastqc/fastqc_data.txt'
    picard = sample + '.picard_RNAseq_qc.txt'
    rd_len = parseFASTQC(fastqc, loc)

    (tot_pairs, r1a_pct, r2a_pct, short, rp_pass, tot_bp, r1_trim, r2_trim, bp_pass) = parseCUTADAPT(cutadapt, loc)

    date_aligned = time.strftime("%c")
    (start_reads, pct_uniq_map, annot_sj, new_sj, pct_mm_pb, pct_del_bp, pct_ins_pb, pct_mm, pct_mmm, pct_unmapped)\
        = parseSTAR(star, loc)
    (median_insert_size, median_absolute_deviation, mean_insert_size, insert_standard_deviation) = parseINS(insert, loc)
    picard_rnaseq_dict = parsePICARD(picard, loc)

    json_dict = {'BionimbusID': RG[0], 'Date': RG[1], 'Machine': RG[2], 'Run': RG[3], 'BarCode': RG[4],
                 'Lane': RG[5], 'read_length': rd_len, 'strand': strand, 'align_date': date_aligned,
                 'cutadapt_stats': {'r1_adapter': r1, 'r2_adapter': r2, 'min_len': minlen, 'min_qual': minqual,
                                    'starting_read_pairs': tot_pairs, 'pct_r1_adapt': r1a_pct, 'pct_r2_adapt': r2a_pct,
                                    'pct_too_short': short, 'rp_pass': rp_pass, 'total_bp': tot_bp,
                                    'pct_r1_qtrim': r1_trim, 'pct_r2_qtrim': r2_trim, 'pct_bp_passed': bp_pass},
                 'STAR_stats': {'input_reads': start_reads, 'pct_uniq_map': pct_uniq_map, 'annot_sj': annot_sj,
                                'non-canon_sj': new_sj, 'mismatch_per-base': pct_mm_pb, 'del_per-base': pct_del_bp,
                                'ins_per-base': pct_ins_pb, 'pct_multi-map': pct_mm, 'pct_uber-multi-map': pct_mmm,
                                'pct_unmapped': pct_unmapped},
                 'picard_stats': {'med_insert_size': median_insert_size, 'med_abs_dev': median_absolute_deviation,
                                  'x_ins_size': mean_insert_size, 's_ins_size': insert_standard_deviation}}
    json_dict['picard_stats'].update(picard_rnaseq_dict)

    json_out = open(sample + '.qc_stats.json', 'w')
    json_out.write(json.dumps(json_dict, indent=4, sort_keys=True))
    json_out.close()

    return 0

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Coverage and algnment QC summary.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sa', '--sample', action='store', dest='sample', help='Sample prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample) = (inputs.config_file, inputs.sample,)
    parse_qc(config_file, sample)