RNAseq Paired End pipeline
===========================
Adapted from Xiyao's pipeline to run on PDC by Miguel Brown, 2015 May

## UTILITY:

#### basic_node_setup.py
usage: basic_node_setup.py [-h] [-id BID] [-j CONFIG_FILE] [-w WAIT]

VM spawner for pipeline processes. Sets up vm for sample analysis and attaches
storage references

optional arguments:
  -h, --help            show this help message and exit
  -id BID, --BID BID    Project Bionimbus ID
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file with snapshot ids and set up params
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Reommended value 600 (in seconds)

#### date_time.py
Simple helper module that prints the current timestamp
#### attach_cinder.py
usage: attach_cinder.py [-h] [-sid SID] [-vid VID] [-id BID] [-s SIZE]
                        [-ip IP] [-w WAIT]

Attaches cinder volume with references to existing vm

optional arguments:
  -h, --help            show this help message and exit
  -sid SID, --snapshot-id SID
                        ID of snapshot. Use cinder to find
  -vid VID, --virt-mach VID
                        Virtual machine id. Use Nova to find
  -id BID, --BID BID    Bionimbpus project id
  -s SIZE, --size SIZE  Cinder reference size. Recommended value 200 (in GB)
  -ip IP, --ip_add IP   VM IP address
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Recommended value 300 (in seconds)

#### cleanup.py                                                                                                                                               
usage: cleanup.py [-h] [-cid CID] [-vid VID] [-id BID] [-ip VIP]

Breaks down a vm built with the standard of having the bionimbus project ID as
part of it's name and attached to a reference

optional arguments:
  -h, --help            show this help message and exit
  -cid CID, --cinder-id CID
                        ID of attached cinder volume. Use cinder to find
  -vid VID, --virt-mach VID
                        Virtual machine id. Use Nova to find
  -id BID, --BID BID    Bionimbpus project id
  -ip VIP, --ip_add VIP
                        VM IP address

##### config_files/GENCODE19_config.json
JSON config file with standard references and tools locations

#### mount.sh
Command basic_node_setup.py uses to mount a reference to a vm.

#### setup_vm.py
usage: setup_vm.py [-h] [-id BID] [-im IMAGE] [-w WAIT] [-f FLAVOR] [-k KEY]

VM spawner for pipeline processes

optional arguments:
  -h, --help            show this help message and exit
  -id BID, --BID BID    Project Bionimbus ID
  -im IMAGE, --image IMAGE
                        Image id to spawn
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Reommended value 300 (in seconds)
  -f FLAVOR, --flavor FLAVOR
                        Image "flavor" to spawn
  -k KEY, --key KEY     Image key-pair to use

##### std_vm_config.json
JSON configuration parameters for creating a pipeline vm and attaching reference storage to it

#### unmount.sh
Command cleanup.py uses to unmount a reference from a vm.

#### download_from_swift.py 
usage: download_from_swift.py [-h] [-c CONT] [-o OBJ]

Simple download module to get files from swift. Can use prefix or whole object
name

optional arguments:
  -h, --help            show this help message and exit
  -c CONT, --container CONT
                        Swift container, i.e. PANCAN
  -o OBJ, --object OBJ  Swift object name/prefix, i.e. RAW/2015-1234

#### upload_variants_to_swift.py 
usage: upload_variants_to_swift.py [-h] [-o OBJ] [-c CONT] [-sl SAMPLE_LIST]
                                   [-sp SAMPLE_PAIRS]

Uploads current directory contents to specified object and container

optional arguments:
  -h, --help            show this help message and exit
  -o OBJ, --object OBJ  Swift object name root to use for aligned merged bam
                        files. i.e. ALIGN/2015-1234
  -c CONT, --container CONT
                        Swfit container name to upload to. i.e. PANCAN
  -sl SAMPLE_LIST, --sample_list SAMPLE_LIST
                        Sample list, one per line
  -sp SAMPLE_PAIRS, --sample_pairs SAMPLE_PAIRS
                        Sample tumor/normal pairs, tsv file with bid pair,
                        sample1, sample2

#### upload_to_swift.py 
usage: upload_to_swift.py [-h] [-o OBJ] [-c CONT]

Uploads current directory contents to specified object and container

optional arguments:
  -h, --help            show this help message and exit
  -o OBJ, --object OBJ  Swift object name to upload current directory contents
                        to. i.e. ALIGN/2015-1234
  -c CONT, --container CONT
                        Swfit container to upload to. i.e. PANCAN

#### update_couchdb.py 
usage: update_couchdb.py [-h] [-f FN]

Update couch db with qc stats using a json object list file

optional arguments:
  -h, --help        show this help message and exit
  -f FN, --file FN  qc_stats.json document list

#### synapse_upload.py
#####Script to upload data for PSYCHencode project.  Could be repurposed or used as as template.
Usage:  ./synapse_upload.py {fixed fields list}{variable fields table}{container}

## ALIGNMENT:
#### pipeline.py 
usage: pipeline.py [-h] [-f1 END1] [-f2 END2] [-t SEQTYPE] [-j CONFIG_FILE]
                   [-m REF_MNT]

RNA alignment paired-end QC pipeline

optional arguments:
  -h, --help            show this help message and exit
  -f1 END1, --file1 END1
                        First fastq file
  -f2 END2, --file2 END2
                        Second fastq file
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file containing tool and reference
                        locations
  -m REF_MNT, --mount REF_MNT
                        Drive mount location. Example would be
                        /mnt/cinder/REFS_XXX
##### Runs the following submodules in order:
1. cutadapter
2. bwt2_pe
3. novosort_sort_pe
4. fastqc
5. picard_insert_size 
6. star


## Pipeline submodule descriptions:

#### cutadapter.py [-h] [-sa SAMPLE] [-f1 END1] [-f2 END2] [-j CONFIG_FILE]

cutadapt module. Removes 3' adapters and trims bases if necessary.Also can
enforce minimum read length - 15 recommended

optional arguments:
  -h, --help            show this help message and exit
  -sa SAMPLE, --sample SAMPLE
                        Sample/location name prefix
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file containing tool and reference
                        locations

#### bwt2_pe.py 
usage: bwt2_pe.py [-h] [-b BWT_TOOL] [-br BWT_REF] [-f1 END1] [-f2 END2]
                  [-s SAMTOOLS_TOOL] [-sr SAMTOOLS_REF] [-sa SAMPLE] [-t T]
                  [-l LOG_DIR]

Bowtie2 paired-end alignment module. Used for insert size estimation toward
beginning using read subset.

optional arguments:
  -h, --help            show this help message and exit
  -b BWT_TOOL, --bwt BWT_TOOL
                        Location of bowtie2 alignment tool.
  -br BWT_REF, --bwt_reference BWT_REF
                        Location of bwt reference file
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file
  -s SAMTOOLS_TOOL, --samtools SAMTOOLS_TOOL
                        Location of samtools tool. Version 1.19 preferred.
  -sr SAMTOOLS_REF, --samtools_reference SAMTOOLS_REF
                        Location of samtools reference
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -t T, --threads T     Number of threads
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

#### novosort_sort_pe.py 
usage: novosort_sort_pe.py [-h] [-n NOVOSORT] [-sa SAMPLE] [-l LOG_DIR] [-t T]
                           [-m MEM]

novosort tool to sort BAM module.

optional arguments:
  -h, --help            show this help message and exit
  -n NOVOSORT, --novosort NOVOSORT
                        novosort binary location
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location
  -t T, --threads T     Number of threads
  -m MEM, --memory MEM  Memory - in GB

#### fastqc.py
usage: fastqc.py [-h] [-f FASTQC_TOOL] [-sa SAMPLE] [-f1 END1] [-f2 END2]

fastqc module. Provides quality stats for fastq file and is independent of
alignment.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTQC_TOOL, --fastqc FASTQC_TOOL
                        Location of fastqc tool.
  -sa SAMPLE, --sample SAMPLE
                        Sample/location name prefix
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file

#### picard_insert_size.py
usage: picard_insert_size.py [-h] [-j JAVA_TOOL] [-p PICARD_TOOL] [-sa SAMPLE]
                             [-l LOG_DIR]

Picard collect insert size metrics module. Gathers insert size metrics, run
after removing BAM duplicates.

optional arguments:
  -h, --help            show this help message and exit
  -j JAVA_TOOL, --java JAVA_TOOL
                        Java location directory, version jdk1.7.0_45 preferred
  -p PICARD_TOOL, --picard PICARD_TOOL
                        Picard jar file location
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

####star.py 
usage: star.py [-h] [-s STAR] [-g GENOME] [-f1 END1] [-f2 END2] [-sa SAMPLE]
               [-l LOG_DIR] [-th TH]

star paired-end alignment module.

optional arguments:
  -h, --help            show this help message and exit
  -s STAR, --star STAR  Location of star aligner.
  -g GENOME, --genome GENOME
                        Location of directory containing genome built by STAR
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location
  -th TH, --threads TH  Number of threads

### Scripts from tophat-cufflinks pipeline - can be used, but are not run in this pipeline script
#### tophat.py
usage: tophat.py [-h] [-t TOPHAT_TOOL] [-tx TX] [-b BWT2_REF] [-f1 END1]
                 [-f2 END2] [-x X] [-sd S] [-sa SAMPLE] [-l LOG_DIR] [-th TH]

tophat paired-end alignment and transcript assembly module.

optional arguments:
  -h, --help            show this help message and exit
  -t TOPHAT_TOOL, --tophat TOPHAT_TOOL
                        Location of tophat alignment tool.
  -tx TX, --transcriptome TX
                        Location of pre-built transcriptome
  -b BWT2_REF, --bwt2_reference BWT2_REF
                        Location of bowtie2 reference file
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file
  -x X, --mean X        Mean insert size
  -sd S, --standard_deviation S
                        Standard deviation of insert size
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location
  -th TH, --threads TH  Number of threads

#### align_stats.py
usage: align_stats.py [-h] [-sa SAMPLE]

Alignment summary report. Converts tophat alignment summary to a table format

optional arguments:
  -h, --help            show this help message and exit
  -sa SAMPLE, --sample SAMPLE
                        Sample/location name prefix

#### cufflinks.py 
usage: cufflinks.py [-h] [-c CUFFLINKS_TOOL] [-e ENS_REF] [-g GENOME]
                    [-sa SAMPLE] [-l LOG_DIR] [-t T]

Annotation and fpkm estimation. Run after tophat.

optional arguments:
  -h, --help            show this help message and exit
  -c CUFFLINKS_TOOL, --cufflinks CUFFLINKS_TOOL
                        Location of cufflinks tool.
  -e ENS_REF, --ensembl_reference ENS_REF
                        Location of ensembl reference file
  -g GENOME, --genome GENOME
                        Location of genome reference file
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location
  -t T, --threads T     Number of threads

#### report.py 
usage: report.py [-h] [-sa SAMPLE] [-r REF_GTF] [-c TX_GTF]

Transcript annotation summart report. Parses cufflinks trnascripts.gtf and
adds gene symbols and transcript biotype

optional arguments:
  -h, --help            show this help message and exit
  -sa SAMPLE, --sample SAMPLE
                        Sample/location name prefix
  -r REF_GTF, --reference REF_GTF
                        Reference gtf file used to annotate - preferably
                        GENCODE with gene name and genetype
  -c TX_GTF, --cufflinks TX_GTF
                        Cufflinks transcript output file

## ANALYSIS
#### ./merge_tables.py -h
Merges tables based on a single column, column denoted array-style. Designed for variant output - will merge gene, position, and base change as name for row

Usage: ./merge_tables.py (<list> <col> <hflag> <tflag>) [<suffix>]

Arguments:
<list>   list of tables
<col>    index of column to merge, array style
<hflag>  1 if tables have header, 0 if not
<tflag>  1 if showing on-target only
<suffix> if given, will omit from column name.  otherwise file name used as column header

#### ./clust_data.py -h
Clustermap generator.

Usage: 
    clust_data.py (-t TABLE) (-o OUTPUT) [-l LOGNORM] [-z ZERO] [-c CMAP] [-d DESCRIPTIVE]

Options:
    -h, --help
    -t TABLE         table name
    -o OUTPUT        output file name
    -l LOGNORM       set id data is to be log transformed
    -z ZERO          zero value to use when log transforming
    -c CMAP          hex colormap to import, otherwise default used
    -d DESCRIPTIVE   add second table with descriptive data to create
                     seperate heatmaps clustered based on data. Ids should be in rows, with header for each row in first column 


## ANNOTATION
#####Rough scripts to add gene names and RNA type to read counting and transcript quantification scripts.  

####annot_htseq.py

####annot_star.py

####gtf2index.py

####kout2ann.py