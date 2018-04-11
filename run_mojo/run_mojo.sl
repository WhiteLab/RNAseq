#!/usr/bin/env bash

export LD_LIBRARY_PATH=/cephfs/users/mbrown/PIPELINES/TOOLS/MOJO/MOJO-v0.0.5-linux-x86_64/lib\:$LD_LIBRARY_PATH
export PATH=/cephfs/users/mbrown/PIPELINES/TOOLS/MOJO/MOJO-v0.0.5-linux-x86_64/bin\:$PATH

$mojo --sample $bnid --json $j --fastq1 $a --fastq2 $b