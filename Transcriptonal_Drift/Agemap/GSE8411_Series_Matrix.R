# These commands take an AGEMAP GSEMatrix File and turn it into an ExpressionSet Object
# GSE8411 file contains z-scores for the Hippocampus
library(GEOquery)
gsm <- getGEO(filename="GSE8411_series_matrix.txt.gz", GSEMatrix=TRUE)
gse8411_expression_matrix <- exprs(gsm)
# Call functions for description of the AGEMAP GSEMatrix file
show(gsm)
dim(gse8411_expression_matrix)
colnames(gse8411_expression_matrix)
length(rownames(gse8411_expression_matrix))
# These commands will extract phenotype data out of AGEMAP GSEMatrix files from pData object
dim(pData(gsm))
rownames(pData(gsm))
colnames(pData(gsm))
# The colnames and row names will vary, so these steps below are needed to know which columns to extract
# pData(gsm)[1,1:10]
# pData(gsm)[1,6]
# pData(gsm)[[1,2]]
# pData(gsm)[[1,1]]
# pData(gsm)[[2,1]]
# sample_title contains tissue, array id, and sample title
sample_title <- as.list(levels(pData(gsm)[[1,1]]))
# sample_id contains the geo_accession numbers for each individual array sample
sample_id <- as.list(levels(pData(gsm)[[1,2]]))
# sample_description contains sample title information as a list (will use this to ensure the order of the samples and indentifiers does not change)
sample_description <- as.list(levels(pData(gsm)[[1,20]]))
# There was a duplicate for first two samples, GSM208884 info missing but was present in sample_id list
# Likely data input error, not a duplicate sample, so corrected sample_description
second_twentyseventh <- sample_description[2:27]
twentyeigth_sixtyseventh <- sample_description[28:67]
twentyseventh_GSM208884 <- "27 AMC76 M 714 A"
sample_descriptor_list <- c(second_twentyseventh,twentyseventh_GSM208884,twentyeigth_sixtyseventh)
# Convert list of character strings to vector of character strings
sample_descriptor_vector <- unlist(sample_descriptor_list)
# Convert vector of character strings into data frame with columns corresponding to each element of string seperated by whitespace
sample_metadata <- read.table(text = sample_descriptor_vector, sep = " ", colClasses = "character")
colnames(sample_metadata) <- c("sample_no","animal_id","sex","age_days","diet")
# Converting sample_id to a vector from a list of strings that contains geo_accession numbers for each sample, and putting that data as a column in the sample_metadata data frame
sample_metadata$geo_accession <- unlist(sample_id)
# Converting sample_title to vector from a list of strings, converting to data frame with columns corresponding to each element of string seperated by comma, and merging the data frame to the sample_metadata data frame
sample_title_vector <- unlist(sample_title)
sample_title_df <- read.table(text = sample_title_vector, sep = ",", colClasses = "character")
colnames(sample_title_df) <- c("tissue","array_id","sample_title")
phenoData <- cbind(sample_title_df, sample_metadata)
age_days <- as.integer(phenoData$age)
age_weeks <- age_days/7
age_months <- age_weeks/4
phenoData$age_months <- age_months
males_adlib <- subset(phenoData,sex=="M" & diet=="A", select=c(1:10))
males_adlib_sample_no <- as.integer(males_adlib$sample_no)
females_adlib <- subset(phenoData,sex=="F" & diet=="A", select=c(1:10))
females_adlib_sample_no <- as.integer(females_adlib$sample_no)
adlib_group <- subset(phenoData, diet=="A", select=c(1:10))
adlib_sample_no <- as.integer(adlib_group$sample_no)
males_cr <- subset(phenoData,sex=="M" & diet=="C", select=c(1:10))
males_cr_sample_no <- as.integer(males_cr$sample_no)
females_cr <- subset(phenoData,sex=="F" & diet=="C", select=c(1:10))
females_cr_sample_no <- as.integer(females_cr$sample_no)
cr_group <- subset(phenoData, diet=="A", select=c(1:10))
cr_group_sample_no <- as.integer(cr_group$sample_no)