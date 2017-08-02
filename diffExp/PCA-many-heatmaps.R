# allows for use with Rscript bash functionality
args <- commandArgs(TRUE)
# if using half-rack VM, enviromental proxy setting is required
Sys.setenv(http_proxy="http://cloud-proxy:3128")
Sys.setenv(https_proxy="http://cloud-proxy:3128")
options(BioC_mirror="http://master.bioconductor.org")
source("http://master.bioconductor.org/biocLite.R")

#prcomp is used instread of princomp here to perform PCA since:
# 1) princomp is limited to experiments where observations >> variables
# 2) prcomp can be used to easily visualize clusters using PCs and metadata

# ggfortify and ggplot are required for plotting and clustering PCs to metadatas
#library(ggfortify)
library(ggplot2)
library(data.table)
library(gplots)

# output_file is the name of the PDF the user would like to create and write to
# the pdf() function actually creates the pdf
# temp args array to test script
#args = c("U049MAI_pca.pdf", "U049MAI_gene_eff_cts.csv", "test_metadata.csv", "cat_list.txt", "cont_list.txt", 6)
output_file = args[1]
pdf(file = output_file)

# counts_file is the FULL path to the counts matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             counts must be rounded to nearest integer
#             row names should be gene/transcript names
#             headers should be BIDs or sample ID names
# meta_file is the FULL path to the metadata matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             row names should be BIDs or sample ID names
#             headers should be the name of the metadata measurement (i.e. sex, diagnosis, chrM contamination, etc...)
# to keep it straight, it makes sense to have n x m counts files and a m x u because the matrices can be multiplied since
# their inner components are the same dimensions ( nxm matrix [dotproduct] mxu matrix = nxu matrix)
counts_file = args[2]
meta_file = args[3]
# give two variable list files pertaining to which are categorical and which are continuous
cat_list <- scan(args[4], what="", sep="\n")
cont_list <- scan(args[5], what="", sep="\n")
# set stopping point for
sub_stop = args[6]

# read in gene expression table and metadata row.names needs to be set to the column number where the row names are listed
# it is important to set check.names=FALSE so that R doesn't try to change any numerical names into odd characters
total_raw_counts <- read.table(counts_file, header=TRUE, row.names=1, check.names=FALSE, sep=",")
metadata <- read.table(meta_file, header=TRUE, row.names=1, check.names=FALSE, sep=",")

# converts each table into a data frame
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

#pre-filter step, remove any ROWS that have zero counts in the count matrix, it does not mean that a sample cannot
# have zero counts, we are just removing the genes where no counts exist across all samples
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe[rowSums(total_raw_read_counts_dataframe)!=0, ]

#add pseudo-count, why??  Because the log(0) = Does not exist, so to address this issue add 1 to all counts
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe + 1
# take the log of the counts data, this will help normalize the data
transformed_dataframe <- log(total_raw_read_counts_dataframe)

# remove BIDs that have missing or strange data
# transformed_dataframe <- within(transformed_dataframe, rm('nothing'))
# metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('nothing'), ]

# The remaining block of code can be uncommented out if you are going to use a SUBSET of the data for PCA
# in the metadata_to_be_removed variable you can indicate in which() statment the column with the specified variable
# of all the samples to remove for analysis.  For example, in this case we go to the column 'Diagnosis' and remove
# all samples that are listed as 'Control'.  This can be changed to anything.
# metadata_to_be_removed <- metadata_dataframe[which(metadata_dataframe$Diagnosis=='SCZ'), ]
# sample_removal <- levels(droplevels(metadata_to_be_removed$BID))
# transformed_dataframe <-transformed_dataframe[,-which(names(transformed_dataframe) %in% sample_removal)]
# metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% sample_removal, ]
# print("Samples to be removed from counts matrix and metadata")
# print(sample_removal)


#remove rows that have a sum of zero
transformed_dataframe <- transformed_dataframe[rowSums(transformed_dataframe)!=0, ]

# sort dataframes.  Dataframes MUST be sorted by colnames in the counts matrix
# and sorted by the rownames in the metadata matrix.  This ensures that the sample names properly match
# each other between counts matix and metadata matrix.  Note, you will see the word 'TRUE' printed
# if they are properly match, else you will see 'FALSE' in which case you need to sort
counts_sorted <- transformed_dataframe[,order(colnames(transformed_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]
all(rownames(metadata_sorted)==colnames(counts_sorted))

# This part is more important if you are removing data from analysis.  It has to extract and collapse missing
# levels from the dataframe of data that was removed
print("getting levels")
sapply(metadata_dataframe,levels)

# Calculation of the prinipal components using prcomp
# the t() function means to take the transpose of the counts matrix
# scale.=TRUE means PCs will be baed on corrlation matrix and not covariance matrix
# corr is better if scales of variables are very different
# center=TRUE uses the centers the data, use centroid/mean
pca_matrix <- prcomp(t(counts_sorted), center=TRUE, scale. = TRUE)

# plot the PCs againts how much variance each PC contributes to the data set
# looking to incorporate the min number of PCs before "elbowing-effect" into the model unless
# a PC is strongly correlated with a variable in the metadata set, in which case,
# just regress out the variable rather than the PC
plot(pca_matrix, type ="l")

#********************************************Model functions***************************************
# helper function to iterate through desired fields to regress in two modes - categorical and continuous
build_categorical_model <- function(factors_affecting_pcs, category, metadata_sorted, pca_matrix, pc){
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,category])))
    factors_affecting_pcs[[category]][[as.character(pc)]]=list()
    factors_affecting_pcs[[category]][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[[category]][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1])
    return(factors_affecting_pcs)
}

build_continuous_model <- function(factors_affecting_pcs, continuous_variable, metadata_sorted, pca_matrix, pc){
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,continuous_variable]))
    factors_affecting_pcs[[continuous_variable]][[as.character(pc)]]=list()
    factors_affecting_pcs[[continuous_variable]][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[[continuous_variable]][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1])
    return(factors_affecting_pcs)
}
# ************************End model functions*****************************************************
# ***************************Function to correlate PCs with metadata******************************
model_pcs <- function(pca_matrix){
  # this block of code will need to be changed depending on the metadata columns available
  factors_affecting_pcs=list()
  # iterate through all PCs and get -log10(p-value) and adjusted R-squared value to every PC correlated to each
  # metadata column listed
  for (pc in seq(1,dim(pca_matrix$rotation)[2])){
    # correlation PCs on factor sex
    # lm() is the Limma library function to build a linear model, in this case each PC is being tested against 'Sex'
    # to determine if any PC is strongly correlated to Sex.  If yes, regress sex out in model
    # Note: in the linear model as.factor() should be placed around categorical variables, while it can be omitted
    # in continuous non-discrete variables
    for (category in cat_list){
        factors_affecting_pcs = build_categorical_model(factors_affecting_pcs, category, metadata_sorted, pca_matrix, pc)
    }
    for (continuous_variable in cont_list){
        factors_affecting_pcs = build_continuous_model(factors_affecting_pcs, continuous_variable, metadata_sorted, pca_matrix, pc)
    }
  }

  # create heatmeap to visualize PC and metadata correlations
  # create a dataframe to store all the -log10 p-values and adjusted R squared vals for the visualization of PCA data
  pvalues <- data.frame()
  adjRsq <- data.frame()

  # iterate only through the look first 10 PCs and extract their -log10(pval) and adj R-sq values
  for (all_factors in seq(1,length(factors_affecting_pcs))){
    for (pc in seq(1,sub_stop)){
      pvalues[all_factors, pc] <- unlist(factors_affecting_pcs[all_factors][[1]][[pc]][[2]])
      adjRsq[all_factors, pc] <- unlist(factors_affecting_pcs[all_factors][[1]][[pc]][[1]])
    }
  }

  # get the row and column names match the p-values
  rownames(pvalues) <- names(factors_affecting_pcs)
  colnames(pvalues) <- unlist(lapply(seq(1,sub_stop),function(x) paste(c('PC',x),collapse='')))

  # get the row and column names matching the adj R-sq values
  rownames(adjRsq) <- names(factors_affecting_pcs)
  colnames(adjRsq) <- unlist(lapply(seq(1,sub_stop),function(x) paste(c('PC',x),collapse='')))

  # round all -log10(pvalue) in the dataframe to three decimal places
  is.num <- sapply(pvalues, is.numeric)
  pvalues[is.num] <- lapply(pvalues[is.num], round, 3)

  # create a heatmap of these values, value is -log10(p-val) and color is the adj R-sq value
  heatmap.2(as.matrix(adjRsq), cellnote=pvalues, notecol = "black", notecex = 0.5, cexRow = 0.3, dendrogram = "none", col=colorRampPalette(c("white", "yellow", "red"))(10))
  print("heatmap completed")
}
#*********************************************End of function*******************************************

# ********************************************Regression functions**************************************
regress_categorical <- function(pca_matrix, metadata_sorted, category){
    for (pc in seq(1,dim(pca_matrix$rotation)[2])){
        linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,category])))
        pca_matrix$x[,pc]  <- linear_model$residuals
    }

    # call function again on regressed out variables
    model_pcs(pca_matrix)
    # writes the variables that were regressed out, NOTE this must be changed manually!!!!!!
    mtext(paste0("vars regressed: ", category), side=3, line=0)

}
regress_continuous <- function(pca_matrix, metadata_sorted, continuous){
    for (pc in seq(1,dim(pca_matrix$rotation)[2])){
      linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,continuous]))
      pca_matrix$x[,pc]  <- linear_model$residuals
  }

    # call function again on regressed out variables
    model_pcs(pca_matrix)
    # writes the variables that were regressed out, NOTE this must be changed manually!!!!!!
    mtext(paste0("vars regressed: ", continuous), side=3, line=0)
}

# ********************************************End regression functions**********************************


model_pcs(pca_matrix)

#---------------------------regress out variables of interest-------------------------------------------------
# for each variable that is to be regressed a linear model must be made and the residuals of the linear model
# must be extract and replace the old PC matrix
# regress out FlowcellBatch
for (category in cat_list){
    regress_categorical(pca_matrix, metadata_sorted, category)
}
for (cont in cont_list){
    regress_continuous(pca_matrix, metadata_sorted, cont)
}

# saves and closes newly created PDF
dev.off()

# call function again on regressed out variables
#model_pcs(pca_matrix)
# writes the variables that were regressed out, NOTE this must be changed manually!!!!!!
#mtext("vars regressed: FlowcellBatch, Sex, UF 5-3 bias, RIN, PMI, AgeDeath, TissueState, BrainBank", side=3, line=0)
# saves and closes newly created PDF
#dev.off()
