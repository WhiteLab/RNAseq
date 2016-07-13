
    cur_mean = int(round(mean(means)))
    cur_std = int(round(mean(stds)))
    novosort_merge_pe(inputs.config_file, bid)
    check = express_quant(bid, inputs.config_file, ref_mnt, str(cur_mean), str(cur_std))
    if check != 0:
        log(loc, date_time() + 'Quantification of RNA failed.  Please check logs\n')
        exit(1)
    mv_cmd = 'mv *.bam BAMS/; mkdir REPORTS; mv *xpr* REPORTS/;'
    subprocess.call(mv_cmd, shell=True)
    log(loc, date_time() + 'Uploading merged bam and quant files for ' + bid + '\n')
    upload_special(bid, cont, obj)
    os.chdir('../../')