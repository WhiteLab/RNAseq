#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/RNAseq/bin/activate
$pipeline --file1 $f1 --file2 $f2 --json $j