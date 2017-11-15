#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu

$express $transcriptome $bam --no-update-check -o $wd --$strand -m $x -s $s --logtostderr 2>> $loc
mv $wd'results.xprs' $wd$root'.express_quantification.txt'
mv $wd'params.xprs' $wd$root'.params.xprs'