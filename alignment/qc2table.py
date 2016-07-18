#!/usr/bin/env python

import json
import os
import sys
sys.path.append('/home/ubuntu/TOOLS/Scripts/utility')
from date_time import date_time
import subprocess
from subprocess import check_output
from log import log


def download_from_swift(cont, obj, lane_list):
    src_cmd = ". /home/ubuntu/.novarc;"
    lanes = open(lane_list, 'r')
    # header doubles as keys for qc stats dicts
    gen_keys = ("BionimbusID", "Date", "Machine", "Run", "BarCode", "Lane", "align_date", "read_length", "strand")
    cutadapt_keys = ("r1_adapter", "r2_adapter", "min_len", "min_qual", "starting_read_pairs", "pct_r1_adapt",
                     "pct_r2_adapt", "pct_too_short", "rp_pass", "total_bp", "pct_r1_qtrim", "pct_r2_qtrim",
                     "pct_bp_passed")
    star_keys = ("input_reads", "pct_uniq_map", "pct_multi-map", "pct_uber-multi-map", "pct_unmapped",
                 "mismatch_per-base", "del_per-base", "ins_per-base", "annot_sj", "non-canon_sj")
    picard_keys = ( "x_ins_size", "s_ins_size", "CORRECT_STRAND_READS", "INCORRECT_STRAND_READS", "CODING_BASES",
                      "INTRONIC_BASES", "INTERGENIC_BASES", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS",
                      "MEDIAN_5PRIME_TO_3PRIME_BIAS", "MEDIAN_CV_COVERAGE", "PCT_CORRECT_STRAND_READS", "UTR_BASES",
                      "PCT_MRNA_BASES", "PCT_INTRONIC_BASES", "PCT_CODING_BASES", "PCT_INTERGENIC_BASES",
                      "PCT_RIBOSOMAL_BASES", "RIBOSOMAL_BASES")
    sys.stdout.write('\t'.join((gen_keys + cutadapt_keys + star_keys + picard_keys)) + '\n')
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = obj + '/' + bid + '/QC/' + bid + '_' + lane + '.qc_stats.json'
            swift_cmd = src_cmd + "swift download " + cont + " --skip-identical --prefix " + cur
            sys.stderr.write(date_time() + swift_cmd + "\n")
            try:
                check = check_output(swift_cmd, shell=True, stderr=subprocess.PIPE)
            except:
                sys.stderr.write(date_time() + "Download of " + obj + " from " + cont + " failed\n")
                exit(1)
            qc_dict = json.loads(open(cur, 'r').read())
            data = []
            for entry in gen_keys:
                data.append(qc_dict[entry])
            for entry in cutadapt_keys:
                data.append(qc_dict['cutadapt_stats'][entry])
            for entry in star_keys:
                data.append(qc_dict['STAR_stats'][entry])
            for entry in picard_keys:
                data.append(qc_dict['picard_stats'][entry])
            sys.stdout.write('\t'.join(data) + '\n')
    lanes.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of qc stats')
    parser.add_argument('-c', '--container', action='store', dest='cont', help='Swift container prefix, i.e. PANCAN')
    parser.add_argument('-o', '--object', action='store', dest='obj',
                        help='Swift object name/prefix, i.e. RAW/2015-1234')
    parser.add_argument('-l', '--lane_list', action='store', dest='lane_list',
                        help='Original lane list used to run pipeline')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (cont, obj, lane_list) = (inputs.cont, inputs.obj, inputs.lane_list)
    download_from_swift(cont, obj, lane_list)
