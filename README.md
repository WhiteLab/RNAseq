RNAseq Paired End pipeline
===========================
RNAseq analysis pipeline in two parts, alignment QC and quantification.  Tools featured are cutadapt, STAR, eXpress, picard.
Designed to run on an hpc with slurm job controller.
Disclaimer:  There are other modes and tools in progress outside of what's outlined here.  Use at your own risk!

## Quick Start - (ain't nobody got time for that!)
## CRASH COURSE RUN:
You don't have time to read through what each script does, and you're a BAMF.  However, being familiar with the Linux environment is highly recommended.

##### 1) Create job run files
###### a) Get list of fastq files to process from swift, one file per line, using a new-line separated list of bionimbus ids:
```
cat <bnids_list> | xargs -IFN find <fastq file home directory>/BN -name '*.gz' > fastq_list 
```
###### b) Use this list to create a run file, in this example for custom capture:
```
<RNA pipe home dir>/utility/fq2lane.py -f fastq_list -s <seq type> > lane_list
```
Typical fastq list:

2018-254_180309_NB551089_0045_AHJGNFBGX5_1_1_sequence.txt.gz
2018-254_180309_NB551089_0045_AHJGNFBGX5_1_2_sequence.txt.gz
2016-2199_180302_K00216R_0042_BHNVTKBBXX_1_1_sequence.txt.gz
2016-2199_180302_K00216R_0042_BHNVTKBBXX_1_2_sequence.txt.gz
2016-2209_180302_K00216R_0042_BHNVTKBBXX_1_1_sequence.txt.gz
2016-2209_180302_K00216R_0042_BHNVTKBBXX_1_2_sequence.txt.gz

Resultant lane_list:

2018-254	capture	180309_NB551089_0045_AHJGNFBGX5_1
2016-2199	capture	180302_K00216R_0042_BHNVTKBBXX_1
2016-2209	capture	180302_K00216R_0042_BHNVTKBBXX_1

##### c) Check config file - this file is typically in <RNA pipe hoe dir>/utility/config_files/v2.5_config.json, and can be copied and modified.  Fields that are likely to be adjusted:

    "refs":{
    "project_dir": "/cephfs/PROJECTS/",
    "project": "PANCAN", # typically as a subdir if the project dir
    "align_dir": "ALIGN_RNASEQ", # typically as a subdir of project
    "config": "/cephfs/users/mbrown/RNAseq/utility/config_files/slurm_config.json" # full path to config file being used

    },
    "params":{
	"threads":"8",
	"ram":"36",
	"user": "mbrown", # this and next param to set acls
    "group": "CPCI",
    "r1adapt": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", # adapter to trim from read 1
    "r2adapt": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", # adapter to trim from read2
    "strand": "rf-stranded", # strandedness for QC and future qunatification.  change to N if not
    "stranded": "Y" # related flag to above, change to N if not
    }
See example config in utility/config_files/slurm_config.json
##### 2) Pipeline run - QC:

```
<RNAseq pipe dir>/alignment/pipeline_wrapper.py -f lane_list.txt -j modified_config.json 2> run.log
```

The pipeline will iterate throught the list, Create run directories with the align_dir specified above, trim fastqs, align with start, and qc the fastq and bam files.  Logs track most of the steps.  Multiple qc tables can ba concatenated for convenience after run using:
```
<RNAseq pipe dir>/alignment/qc2table.py -h
usage: qc2table.py [-h] [-d PROJECT_DIR] [-p PROJECT] [-a ALIGN_DIR]
                   [-l LANE_LIST]

Uses pipeline lane list to create a summary table of qc stats

optional arguments:
  -h, --help            show this help message and exit
  -d PROJECT_DIR, --project-dir PROJECT_DIR
                        Project dir, i.e. /cephfs/PROJECTS/
  -p PROJECT, --project PROJECT
                        project name, i.e. PANCAN
  -a ALIGN_DIR, --align-dir ALIGN_DIR
                        Alignment subdirectory, i.e. ALIGN_RNASEQ
  -l LANE_LIST, --lane_list LANE_LIST
                        Original lane list used to run pipeline
```
Resultant files are organized into subdirectories: TRIMMED_FQ, BAMS, LOGS, and QC

##### 3) Pipeline run - Quantification:
This will take the same pipeline run files from alignment and run eXpress for transcript-level quantification

```
<RNAseq pipe dir>/analysis/quant_pipe_wrap.py -f lane_list.txt -j modified_config.json 2> run.log
```
It is recommended that you lower the cpu and memory usage as it is less intensive than the alignment pipeline to about half.
eXpress output can be collapsed to a gene-level table, choosing RNA biotype and count type, if desired, using the following helper script:

```
/cephfs/users/mbrown/PIPELINES/RNAseq/utility/collapse_merge_express.py -h
Combines express output tables, collapsing on gene symbol by effective count

Usage: collapse_merge_express.py <ct_list> <type_list> <field> <round> <c_flag> <mode>

Arguments:
  <ct_list> list of express count files
  <type_list> list of transcript types to accept. Type None not to use one
  <field> name of field to collapse on or default for est_counts
  <round> 'y' to round values after collapsing
  <c_flag> 'y' to collapse by gene or just merge by transcript and ENSEMBL id
  <mode> for collapse on gene mode, use 'median' or 'sum'

Options:
  -h 
 ```
#### 4) Optional - fusion detection using MOJO
Ought to be run after alignment QC regardless of the fact that this pipeline uses a completely different set of tools,
details at author's github here: https://github.com/cband/MOJO

```
/cephfs/users/mbrown/PIPELINES/RNAseq/run_mojo/mojo_wrap.py 
usage: mojo_wrap.py [-h] [-l LANE] [-j CONFIG_FILE]

Search for fusion junctions using MOJO

optional arguments:
  -h, --help            show this help message and exit
  -l LANE, --lane LANE  Lane list file
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file containing tool and reference
                        locations

```
See example config in utility/config_files/mojo_wrap_config.json
## Current setup:

#### cutadapt v1.15
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

#### STAR v2.5.2a
Purpose: RNAseq aligner
Download pre-built binary from https://github.com/alexdobin/STAR

#### express 1.5.1
Purpose: Transcript-level quantification
Download from https://pachterlab.github.io/eXpress/
