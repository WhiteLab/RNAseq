#!/usr/bin/env python3

import json
import sys
import os
sys.path.append('/cephfs/users/mbrown/PIPELINES/RNAseq/')


def qc2table(project_dir, project, align_dir, lane_list):
    lanes = open(lane_list, 'r')
    # header doubles as keys for qc stats dicts
    gen_keys = ("BionimbusID", "Date", "Machine", "Run", "BarCode", "Lane", "align_date", "read_length", "strand")
    cutadapt_keys = ("r1_adapter", "r2_adapter", "min_len", "min_qual", "starting_read_pairs", "pct_r1_adapt",
                     "pct_r2_adapt", "pct_too_short", "rp_pass", "total_bp", "pct_r1_qtrim", "pct_r2_qtrim",
                     "pct_bp_passed")
    star_keys = ("input_reads", "pct_uniq_map", "pct_multi-map", "pct_uber-multi-map", "pct_unmapped",
                 "mismatch_per-base", "del_per-base", "ins_per-base", "annot_sj", "non-canon_sj")
    picard_keys = ("x_ins_size", "s_ins_size", "CORRECT_STRAND_READS", "INCORRECT_STRAND_READS", "CODING_BASES",
                   "INTRONIC_BASES", "INTERGENIC_BASES", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS",
                   "MEDIAN_5PRIME_TO_3PRIME_BIAS", "MEDIAN_CV_COVERAGE", "PCT_CORRECT_STRAND_READS", "UTR_BASES",
                   "RIBOSOMAL_BASES", "PCT_MRNA_BASES", "PCT_INTRONIC_BASES", "PCT_CODING_BASES",
                   "PCT_INTERGENIC_BASES", "PCT_RIBOSOMAL_BASES")
    sys.stdout.write('\t'.join((gen_keys + cutadapt_keys + star_keys + picard_keys)) + '\n')
    for line in lanes:
        line = line.rstrip('\n')
        (bid, seqtype, lane_csv) = line.split('\t')
        for lane in lane_csv.split(', '):
            cur = project_dir + project + '/' + align_dir + '/' + bid + '/QC/' + bid + '_' + lane + '.qc_stats.json'
            if os.path.isfile(cur):
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
            else:
                sys.stderr.write('Could not find ' + cur + ', SKIP!\n')
    lanes.close()
    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses pipeline lane list to create a summary table of qc stats')
    parser.add_argument('-d', '--project-dir', action='store', dest='project_dir',
                        help='Project dir, i.e. /cephfs/PROJECTS/')
    parser.add_argument('-p', '--project', action='store', dest='project',
                        help='project name, i.e. PANCAN')
    parser.add_argument('-a', '--align-dir', action='store', dest='align_dir',
                        help='Alignment subdirectory, i.e. ALIGN_RNASEQ')
    parser.add_argument('-l', '--lane_list', action='store', dest='lane_list',
                        help='Original lane list used to run pipeline')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (project_dir, project, align_dir, lane_list) = (inputs.project_dir, inputs.project, inputs.align_dir,
                                                    inputs.lane_list)
    qc2table(project_dir, project, align_dir, lane_list)
