RNAseq Paired End pipeline
===========================
RNAseq analysis pipeline in two parts, alignment QC and qunatification.  Tools featured are cutadapt, STAR, eXpress, picard.
Disclaimer:  There are other modes and tools in progress outside of what's outlined here.  Use at your own risk!

## Quick Start - (ain't nobody got time for that!)
## CRASH COURSE RUN:
You don't have time to read through what each script does, and you're a BAMF.  However, being familiar with OpenStack 
is highly recommended.

##### 1) Set up vm - on head node, will give updates on set up
VMs related to this pipeline have an image with the suffix RNAseqvx.xx and a date prefix.  Best to choose the latest.
Image list can be viewed using:
```
nova image-list
```
Find the ID of the image and boot a vm from it.  You can get a list of "flavors;" this sets up number of cpus, ram, and
 disk space.  An flavor with 8 cpus, 32GB RAM, and 400GB ephemeral is recommended.  Flavors can be listed using the 
 following:
```
nova flavor-list
```
Boot vm with your open stack key:
```
nova boot --image <image-id> --flavor <flavor-id> --key-name your_key DESIRED_VM_NAME
```

After issuing the command, you'll get an ID for your vm.  You can use this to track boot progress, listing the 
current servers:
```
nova list
```
Be sure to transfer your .novarc credential files upon boot.

##### 2) Get reference files from object store
Most can be found in the container MB_TEST, object prefix GENCODE24GRCH37/ should get most of what is needed for alignment.
**Be certain not to download to the root drive.  Ephemeral space is in /mnt, recommended to make a working directory
there**
```
sudo mkdir /mnt/WORK; sudo chown ubuntu:ubuntu /mnt/WORK; swift download MB_TEST --prefix GENCODE24GRCH37/; mv GENCODE24GRCH37
REFS;
```

##### 3) Create job run files
##### a) Get list of fastq files to process from swift, one file per line, using a new-line seprated list of bionimbus ids:
```
/home/ubuntu/utility/bid_swift_list.py -c <swift container> -o <object prefix> -l <bnids list> > fastq_list 
```
##### b) Use this list to create a run file, in this example for custom capture:
```
/home/ubuntu/utility/fq2lane.py -f <fastq_list> -s <seq type> > lane_list
```
Typical fastq list:

RAW/2014-2230/2014-2230_141212_SN1070_0312_BHB5BNADXX_2_1_sequence.txt.gz
RAW/2014-2230/2014-2230_141212_SN1070_0312_BHB5BNADXX_2_2_sequence.txt.gz
RAW/2014-2231/2014-2231_141212_SN1070_0312_BHB5BNADXX_2_1_sequence.txt.gz
RAW/2014-2231/2014-2231_141212_SN1070_0312_BHB5BNADXX_2_2_sequence.txt.gz
RAW/2014-2232/2014-2232_141212_SN1070_0312_BHB5BNADXX_2_1_sequence.txt.gz
RAW/2014-2232/2014-2232_141212_SN1070_0312_BHB5BNADXX_2_2_sequence.txt.gz
RAW/2014-2232/2014-2232_150501_SN1070_0375_AH3L5KBCXX_1_1_sequence.txt.gz
RAW/2014-2232/2014-2232_150501_SN1070_0375_AH3L5KBCXX_1_2_sequence.txt.gz
RAW/2014-2233/2014-2233_141212_SN1070_0312_BHB5BNADXX_2_1_sequence.txt.gz
RAW/2014-2233/2014-2233_141212_SN1070_0312_BHB5BNADXX_2_2_sequence.txt.gz

Resultant lane_list:

2014-2232	polyA	141212_SN1070_0312_BHB5BNADXX_2, 150501_SN1070_0375_AH3L5KBCXX_1
2014-2233	polyA	141212_SN1070_0312_BHB5BNADXX_2
2014-2230	polyA	141212_SN1070_0312_BHB5BNADXX_2
2014-2231	polyA	141212_SN1070_0312_BHB5BNADXX_2

##### c) Check config file - this file is typically in ~/TOOLS/Scripts/utility/config_files/v2.5_config.json, and can be copied and modified.  Fields that are likely to be adjusted:

    "refs":{
	"cont":"PANCAN", # container
	"obj":"ALIGN", # object prefix

	"config":"/cephfs/users/mbrown/PIPELINES/DNAseq/utility/config_files/complete_config.json" # this file location
    },
    "params":{
	"threads":"8",
	"ram":"30",
    "r1adapt": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", # adapter to trim from read 1
    "r2adapt": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", # adapter to trim from read2
     "strand": "rf-stranded", # strandedness for QC and future qunatification.  change to N if not
    "stranded": "Y" # related flag to above, change to N if not
    }

##### 4) Pipeline run - QC:

```
/cephfs/users/mbrown/PIPELINES/DNAseq/alignment/pipeline_wrapper.py -f lane_list.txt -j modified_config.json -m location_of_volume_mount 2> run.log
```
-m clarification:
 -m REF_MNT, --mount REF_MNT
                        Reference drive mount location. Example would be
                        /mnt/cinder/REFS_XXX

The pipeline will iterate throught the list upload files to swift, and delete on the volume for next run.  Logs track most of the steps.  Multiple qc tables can ba concatenated for convenience after run using /cephfs/users/mbrown/PIPELINES/DNAseq/alignment/merge_qc_stats.py:
```
/cephfs/users/mbrown/PIPELINES/DNAseq/alignment/qc2table.py -f <lane_list> -c <swift container> -o <swift object prefix> > qc_table.txt 2> log.txt
```

##### 5) Pipline run - Quantification:
This will take the same pipeline run files from alignment and run eXpress for transcript-level quantification

```
/cephfs/users/mbrown/PIPELINES/DNAseq/analysis/quant_pipe.py -f lane_list.txt -j modified_config.json -m location_of_volume_mount 2> run.log
```

eXpress output can be collapsed to a gene-level table, choosing RNA biotype and count type, if desired, using the following helper script:

```
/cephfs/users/mbrown/PIPELINES/DNAseq/utility/collapse_merge_express_by_gene.py <ct_list> <type_list> <field> > output.txt
```
Arguments:
  <ct_list> list of express count files
  <type_list> list of transcript types to accept. Type None not to use one
  <field> name of field to collapse on or default for est_counts


# Software requirements
If you are starting from scratch, with some modification, this pipeline can be run on a blank vm running linux or a server that does'nt use an object store.  Be aware that new versions of all likely exist now.  You can use them but you'd have to validate the options and functionality in this framework yourself!

## Current setup:

#### cutadapt v1.8.1
Purpose: Adapter and base quality trimming
can be installed using:

```
pip install cutadapt
```

#### bowtie2 2.2.3
Purpose:  Used quickly at the start for insert size estimation
Download and build from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/

#### picard v2.0.1
Purpose: Insert size estimation and RNAseq QC metrics
Download jar from https://broadinstitute.github.io/picard/.  Requires java jdk

#### novosort 1.03.09
Purpose:  Merge and sort bams
Note, this is proprietary software.  Definitely worth it to shave time off of this step as the software is capable of using specified memory and processors to the full extent,  Otherwise, us picard.
Obtain from http://www.novocraft.com/products/novosort/

#### fastqc 0.11.5
Purpose: Fastq QC metrics
Download from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#### STAR v2.4.2a
Purpose: RNAseq aligner
Download pre-built binary from https://github.com/alexdobin/STAR

#### express 1.5.1
Purpose: Transcript-level quantification
Download from https://pachterlab.github.io/eXpress/
