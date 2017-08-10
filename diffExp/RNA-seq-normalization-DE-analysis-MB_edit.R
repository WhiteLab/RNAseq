# allows for use with Rscript bash functionality
args <- commandArgs(TRUE)
if (length(args) < 4){
    stop('Usage: Rscript RNA-seq-normalization-DE-analyis-MB_edit.R {counts.csv} {metadata.csv} {patient_name} {regression values.txt} {control_value}')
}
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
#register(MulticoreParam(7)) # this can be changed to the number of processors/threads to use for parallelization


#-----------------gene expression analysis portion (DE)----------------------

# total_raw_counts MUST have a header and row names and be in CSV format, header=BID or sampleID, row names=gene/transcript names and they must be unique
# metadata MUST have a header and row names and be in CSV format, header=metric names (i.e. sex, diagnosis, RIN, etc...), and row names are the unique sampleIDs or BIDs
# NOTE! row.names=<INTEGER> specifies which column in the data sheets have row names listed, this should be changed depending on the location
# IMPORTANT!  All gene counts MUST BE INTEGERS, no floats, so please round all counts to nearest INT

total_raw_counts <- read.table(args[1], header=TRUE, row.names=1, check.names=FALSE, sep=",")
metadata <- read.table(args[2], header=TRUE, row.names=1, check.names=FALSE, sep=",")
patient = args[3]
# factor_list <- scan(args[4], what="", sep="\n")
control = args[4]
# converts the gene counts matrix and metadata matrix into dataframes
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)


# sorts dataframes so columns and rows are in same BID order for count data and metadata
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# below statement should result in TRUE, else the two dataframes are not properly sorted
all(rownames(metadata_sorted)==colnames(counts_sorted))

# creates a special DESeq object to store the counts data, metadata, and importantly the linear model
# About the linear model:
#   1)  It must always follow the format design=~
#   2)  The condition to be tested MUST ALWAYS BE THE LAST VARIABLE IN THE MODEL, SUPER IMPORTANT!!!!!!!!
#   3)  any covariates in the model that are CATEGORIGAL must have the as.factor() function applied to it or
#       it will treat it as a coninuous variable
#   4)  Any variable that appears in the design formula must be present as a column in the metadata

# hard code factors for now
deseq_obj <- DESeqDataSetFromMatrix(countData= counts_sorted, colData=metadata_sorted , design= ~ as.factor(Kit) + Biotype)

#set comparison reference/control
# remember to change the relevel when comparing non-controls
# this tells DESeq what the baseline sample should be, in this case it the Control labeled patients located in the
# Diagnosis column of the deseq_object
print("Releveling reference to Control")
deseq_obj$Biotype <- relevel(factor(deseq_obj$Biotype), ref=control)
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
# output Cook Distance for outlier detection on a sample level
par(mar=c(8,5,2,2))
boxplot(log10(assays(deseq_obj)[["cooks"]]), range=0, las=2)

# output normalized count data
# write.csv(file="") Input the FULL Path and file name of the normalized counts matrix
norm_counts <- counts(deseq_obj, normalized=TRUE)
write.csv(norm_counts, file=paste(patient, '_norm_counts.csv'))

for (btype in metadata_sorted$Biotype){
    if (btype != control){

        # At this step the results will be extracted.  Parallel refers to the number of processor that will be used to
        # speed up the extraction, if TRUE it means to use parallelization. Contrast tells DESeq how to present the results
        # In this case, it means go to the Diagnosis column of the metadata and put "BP" as the numerator, and "Control"
        # as the denominator when reporting up- and down-regulation of expression.  The essentially means to use "Control"
        # as the baseline level of counts


        de_results_control_bp <- results(deseq_obj, parallel=FALSE, contrast=c("Biotype", btype, control))

        #de_results_control_bp <- results(deseq_obj, parallel=TRUE)

        # This means reorder the results extracted above by padj (FDR) and print out an overall summary to the terminal
        # write.csv(file='') means to write the FULL path and name of the final output file for results
        # A note on what NAs mean in results of csv file:
        # Any gene with NA as p-value and adjusted p-value (FDR) means it is an outlier by Cooks Distance standards
        # in at least 3 or more biological replicates
        resOrdered_c_bp <- de_results_control_bp[order(de_results_control_bp$padj),]
        summary(de_results_control_bp)
        write.csv(as.data.frame(resOrdered_c_bp), file=paste(patient, btype, "vs", control, "_diff_exp.csv", sep="_"))

    }
}