# if using half-rack VM, enviromental proxy setting is required
Sys.setenv(http_proxy="http://cloud-proxy:3128")
Sys.setenv(https_proxy="http://cloud-proxy:3128")
library(data.table)
options(BioC_mirror="http://master.bioconductor.org")
source("http://master.bioconductor.org/biocLite.R")
biocLite(c("DESeq2", "GenomicFeatures", "AnnotationDbi", "limma"))
library("DESeq2")
library("limma")
library("BiocParallel")
register(MulticoreParam(7)) # this can be changed to the number of processors/threads to use for parallelization


#-----------------gene expression analysis portion (DE)----------------------

# total_raw_counts MUST have a header and row names and be in CSV format, header=BID or sampleID, row names=gene/transcript names and they must be unique
# metadata MUST have a header and row names and be in CSV format, header=metric names (i.e. sex, diagnosis, RIN, etc...), and row names are the unique sampleIDs or BIDs
# NOTE! row.names=<INTEGER> specifies which column in the data sheets have row names listed, this should be changed depending on the location
# IMPORTANT!  All gene counts MUST BE INTEGERS, no floats, so please round all counts to nearest INT
total_raw_counts <- read.table("/mnt/WORK/2017-Aug-9_org_proj_rnaseq/study_eff_counts.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
total_raw_counts <- total_raw_counts[ , -which(names(total_raw_counts) %in% c("nothing"))]
metadata <- read.table("/mnt/WORK/2017-Aug-9_org_proj_rnaseq/metadata_clean.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")

# UNCOMMMENT IF ONLY USING A SUBSET OF DATA!! 
# if using a SUBSET of the entire data please provide a metadata sheet with the removed samples
cleaned_specific_metadata <- read.table("/mnt/de_control_bp/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-MISSING-DATA-FILLED-WITH-AVERAGES-add-PC1-control-bp-subset-only-ERCC-removed-from-counts.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
cleaned_specific_metadata_dataframe <-  as.data.frame(cleaned_specific_metadata)

# converts the gene counts matrix and metadata matrix into dataframes
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

# UNCOMMENT IF ONLY USING A SUBSET OF DATA
# remove samples now that you don't want to use for analysis
# change to what samples should be removed
# the which(metadata$Diagnosis=='SCZ') means to remove the SCZ samples located in the Diagnosis column, change to what pertains to the dataset
# the statement sample_removal <- levels(droplevels(metadata_to_be_removed$BID)) means to remove the samples from above matching to the BID column
#--------------------commented out if using entire set-------------------------------------------------------
metadata_to_be_removed <- metadata_dataframe[which(metadata_dataframe$Diagnosis=='SCZ'), ]
sample_removal <- levels(droplevels(metadata_to_be_removed$BID))
total_raw_read_counts_dataframe <-total_raw_read_counts_dataframe[,-which(names(total_raw_read_counts_dataframe) %in% sample_removal)]
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% sample_removal, ]
metadata_sorted <- cleaned_specific_metadata_dataframe[order(rownames(cleaned_specific_metadata_dataframe)),]
print("Samples to be removed from counts matrix and metadata")
print(sample_removal)
#------------------------end of commented out block----------------------------------------------------------

# sorts dataframes so columns and rows are in same BID order for count data and metadata
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
#metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# below statement should result in TRUE, else the two dataframes are not properly sorted
all(rownames(metadata_sorted)==colnames(counts_sorted))

# creates a special DESeq object to store the counts data, metadata, and importantly the linear model
# About the linear model:
#   1)  It must always follow the format design=~
#   2)  The condition to be tested MUST ALWAYS BE THE LAST VARIABLE IN THE MODEL, SUPER IMPORTANT!!!!!!!!
#   3)  any covariates in the model that are CATEGORIGAL must have the as.factor() function applied to it or
#       it will treat it as a coninuous variable
#   4)  Any variable that appears in the design formula must be present as a column in the metadata
deseq_obj <- DESeqDataSetFromMatrix(countData= counts_sorted, colData=metadata_sorted , design= ~ as.factor(Kit) + Biotype)

#set comparison reference/control
# remember to change the relevel when comparing non-controls
# this tells DESeq what the baseline sample should be, in this case it the Control labeled patients located in the
# Diagnosis column of the deseq_object
print("Releveling reference to Control")
deseq_obj$Biotype <- relevel(factor(deseq_obj$Biotype), ref="Normal")

# pre-filtering, remove rows in count data with 0 or 1 reads, this can be changes to anything you want to filter out
# keep if mind for example if you are looking to remove all genes that have only 1 count in 80% of all samples, the 
# 1 listed will have to encorporate a conitional statement to account for this, otherwise it is translated into
# keep all genes where the sums of counts is at greater than one across all samples
deseq_obj <- deseq_obj[rowSums(counts(deseq_obj)) > 1, ]

# perform DE analysis (size factors, dispersion, negative binomial distribution)
# parallel means the process will be parallelized depening on how many processors were specified for use above
# at this step, all total reads will be normalized, variance and dispersion normalization via negative binomial distribution
# will be performed
deseq_obj <- DESeq(deseq_obj, parallel=FALSE, betaPrior=FALSE)

# At this step the results will be extracted.  Parallel refers to the number of processor that will be used to
# speed up the extraction, if TRUE it means to use parallelization. Contrast tells DESeq how to present the results
# In this case, it means go to the Diagnosis column of the metadata and put "BP" as the numerator, and "Control"
# as the denominator when reporting up- and down-regulation of expression.  The essentially means to use "Control"
# as the baseline level of counts
de_results_control_bp <- results(deseq_obj, parallel=FALSE, contrast=c("Biotype", "Tumor", "Normal"))

#de_results_control_bp <- results(deseq_obj, parallel=TRUE)

# This means reorder the results extracted above by padj (FDR) and print out an overall summary to the terminal
# write.csv(file='') means to write the FULL path and name of the final output file for results
# A note on what NAs mean in results of csv file:
# Any gene with NA as p-value and adjusted p-value (FDR) means it is an outlier by Cooks Distance standards
# in at least 3 or more biological replicates
resOrdered_c_bp <- de_results_control_bp[order(de_results_control_bp$padj),]
summary(de_results_control_bp)
write.csv(as.data.frame(resOrdered_c_bp), file="/mnt/WORK/2017-Aug-9_org_proj_rnaseq/T_v_N_kit_regressed_all.csv")

# output Cook Distance for outlier detection on a sample level
par(mar=c(8,5,2,2))
boxplot(log10(assays(deseq_obj)[["cooks"]]), range=0, las=2)

# output normalized count data
# write.csv(file="") Input the FULL Path and file name of the normalized counts matrix
norm_counts <- counts(deseq_obj, normalized=TRUE)
write.csv(norm_counts, file='/mnt/WORK/2017-Aug-9_org_proj_rnaseq/T_v_N_kit_regressed_norm_cts.csv')
