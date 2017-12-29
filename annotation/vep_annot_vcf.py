#!/usr/bin/env python3
import sys
sys.path.append('/cephfs/users/mbrown/PIPELINES/DNAseq/')
from utility.date_time import date_time
from utility.log import log
import json
import subprocess
from vep_report import gen_report


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['tools']['VEP'], config_data['refs']['vepCache'], config_data['refs']['samtools'],\
           config_data['params']['threads'], config_data['refs']['cadd'], config_data['refs']['vep_cache_version']


def pass_filter(sample, somatic_flag):
    in_fn = sample + '.haplo_somatic.vcf'
    out_fn = sample + '.haplo_somatic.PASS.vcf'
    if somatic_flag == 'N':
        in_fn = sample + '.haplo_filtered.vcf'
        out_fn = sample + '.haplo_filtered.PASS.vcf'
    out = open(out_fn, 'w')
    infile = open(in_fn, 'r')
    for line in infile:
        if line[0] == '#':
            out.write(line)
        else:
            fields = line.split('\t')
            if fields[6] == 'PASS':
                out.write(line)
    infile.close()
    out.close()


def annot_gatk_haplotype(config_file, sample_pairs, ref_mnt, somatic_flag):
    (vep_tool, vep_cache, fasta, threads, cadd, cacheVersion) = parse_config(config_file)
    fasta = ref_mnt + '/' + fasta
    cadd = ref_mnt + '/' + cadd
    vep_cache = ref_mnt + '/' + vep_cache
    # scale back on the forking a bit
    if int(threads) > 2:
        threads = str(int(threads)/2)
    # make flexible for a pairs file
    for pair in open(sample_pairs, 'r'):
        pair = pair.rstrip('\n').split('\t')[0]
        pass_filter(pair, somatic_flag)
        in_vcf = pair + '.haplo_somatic.PASS.vcf'
        out_vcf = pair + '.haplo_somatic.PASS.vep.vcf'
        if somatic_flag == 'N':
            in_vcf = pair + '.haplo_filtered.PASS.vcf'
            out_vcf = pair + '.haplo_filtered.PASS.vep.vcf'

        loc = pair + '.vep.log'
        run_vep = 'perl ' + vep_tool + ' --cache -i ' + in_vcf + ' --vcf -o ' + out_vcf + ' --symbol --vcf_info_field' \
                ' ANN --canonical --html --variant_class --sift both --offline --maf_exac --no_whole_genome --fork ' \
                + threads + ' --fasta ' + fasta + ' --dir_cache ' + vep_cache + ' --plugin CADD,' + cadd \
                  + ' --cache_version ' + cacheVersion + ' 2> ' + loc + ' >> ' + loc
        log(loc, date_time() + run_vep + '\n')
        check = subprocess.call(run_vep, shell=True)
        if check != 0:
            log(loc, date_time() + 'VEP annotation failed for ' + pair + '\n')
            exit(1)
        check = gen_report(out_vcf, pair)
        if check != 0:
            log(loc, date_time() + 'VEP as table failed for ' + pair + '\n')
            exit(1)

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Mutation annotation using variant effect predictor.')
    parser.add_argument('-j', '--json', action='store', dest='config_file',
                        help='JSON config file with tool and reference locations')
    parser.add_argument('-sp', '--sample_pairs', action='store', dest='sample_pairs', help='Sample tumor/normal pairs,'
                                                                                           ' or sample list')
    parser.add_argument('-r', '--ref_mnt', action='store', dest='ref_mnt',
                        help='Reference mount directory, i.e. /mnt/cinder/REFS_XXX')
    parser.add_argument('-sf', '--somatic_flag', action='store', dest='somatic_flag',
                        help='Flag whether or not to deal with somatic or single sample calls')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    inputs = parser.parse_args()
    (config_file, sample_pairs, ref_mnt, somatic_flag) = (
        inputs.config_file, inputs.sample_pairs, inputs.ref_mnt, inputs.somatic_flag)
    annot_gatk_haplotype(config_file, sample_pairs, ref_mnt, somatic_flag)
