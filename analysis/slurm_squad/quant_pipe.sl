#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cephfs/software/samtools-1.6/htslib-1.6/htslib/lib
source /cephfs/users/mbrown/.virtualenvs/RNAseq/bin/activate
$quant --sample $bnid --json $j -x $x -s $s