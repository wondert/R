library(GEOquery)
#
#
# Load GSM file with phenotype and expression data
gsm <- getGEO(filename="~/GSE60862/GSE60862_series_matrix.txt.gz", GSEMatrix=TRUE)
######################### exprs(gsm) == expression matrix ##################################### 
# row.names == ID (from fData(gsm)) 
# colnames == SAMPLE ID
############################ fData(gsm) == feature data #######################################
# ID == number representing probe
# GB_LIST == Gene Nomenclature (Alternative Names, IDs Gene/mRNA/Protein)
# gene_assignment == Gene Annotation (Full Annotation)
# mrna_assignment == mRNA Annotation (Full Annotation)
########################### pData(gsm) == phenotype data ######################################
# geo_accession == SAMPLE ID (geo)
# characteristics_ch1.1 == brain region
# characteristics_ch1.3 == study site
# characteristics_ch1.4 == cause of death
# characteristics_ch1.5 == sex
# characteristics_ch1.6 == age at death
# characteristics_ch1.7 == post-mortem interval for tissue collection (hours)
# characteristics_ch1.9 == brain pH
# characteristics_ch1.10 == rin
# title == SAMPLE ID (truncated), brain region
# source_name_ch1 == tissue status, SAMPLE ID, brain region
# data_processing == data processing information
# description == supercedes geo_accesion???
############################### WE HAVE DATA DESCREPANCY ######################################
## length(pData(gsm)$description) != length(unique(pData(gsm)$description))
## length(unique(pData(gsm)$geo_accession)) != length(unique(pData(gsm)$description))
## length(pData(gsm)$description) == length(pData(gsm)$geo_accession)
# Data for both are offset at V2 == GSM1490161
# pData(gsm)$geo_accession[1123] == V1124 == GSM1491283 
# pData(gsm)$description[1123] == V1124 == GSM1135678
# This entry supercedes previous accession: GSM1135615 
# pData(gsm)[1,2];pData(gsm)[1,30]
# j <- pData(gsm)$geo_accession
# j <- levels(j)
# eset <- exprs(gsm)
# names_gsm <- colnames(eset)
# table(names_gsm == j)  yields 1231 TRUE, length both vectors is 1231
# k <- pData(gsm)$description
# k <- levels(k)[2:933]
# k <- gsub('^.{42}', '', k)
# k <- c('',k)
# Note really sure what to do here, these sets dont match up at all
# Will ingore for now.
#
#
#
############################ Subset Analysis -- Hippocampus ###################################
hippocampus <- subset(pData(gsm), characteristics_ch1.1 == 'brain region: hippocampus', select = c('characteristics_ch1.5', 'characteristics_ch1.6'))
hippocampus$characteristics_ch1.5 <- gsub('^.{5}', '', ll$characteristics_ch1.5)
hippocampus$characteristics_ch1.6 <- as.numeric(gsub('^.{25}', '', ll$characteristics_ch1.6))
male_hippocampus <- hippocampus[hippocampus$characteristics_ch1.5 == 'Male',]     # 92 males, age range 16-91, median 55
female_hippocampus <- hippocampus[hippocampus$characteristics_ch1.5 == 'Female',]   # 30 females, age range 20-102, median 61.5
# Need to find for each tissue/sex the index of these samples by each age grouping
# 
#
#
###############################################################################################
###############################################################################################
###############################################################################################
# Here we collect ID and VALUE expression data from SOFT formatted family files.
# Data are RMA quantile normalized expression values that have been summaried by gene from exon level data with Windsored mean
# Download SOFT formatted family file
# In bash shell filter out lines using sed that start with ^ or # or I or $
# sed '/!platform_table_begin/,/!platform_table_end/d;/^#/d;/^\^/d;/^I/d;/^$/d' GSE60862_family.soft.txt > GSE60862_family.soft.cleaned_cp.txt
# NOTE '!' is still present and will be ignored when importing data into R with comment.char = '!'
# per sample data row count = 26493
# No. unique data sets = 1231
# Make sure to add headers c('ID', 'VALUE') to R data.frame
gse_matrix <- "~/GSE60862/GSE60862_family.soft.cleaned_cp.txt"
#
#
#
# Read array data into R. Collect sample/array id's.
gse_gpl5175_norm <- read.table(gse_matrix, header = TRUE, sep = "\t", blank.lines.skip = TRUE, comment.char = "!", fill = TRUE, skip = 0, nrows = -1)
sample_list_numeric <- c(1:1231)
array_id_5175 <- colnames(gse_gpl5175_norm)[1:1231]
colnames(gse_gpl5175_norm) <- c("ID",sample_list_numeric)
#
#
#
gse <- rbind(gse_gpl738_raw,gse_gpl782_raw)
#
#
#
geo_array_id <- c(array_id_738,array_id_782)
sample_title <- as.list(levels(pData(gsm)[[1,1]]))
sample_id <- as.list(levels(pData(gsm)[[1,2]]))
sample_description <- as.list(levels(pData(gsm)[[1,20]]))
second_twentyseventh <- sample_description[2:27]
twentyeigth_sixtyseventh <- sample_description[28:67]
twentyseventh_GSM <- "27 AMC76 M 714 A"
sample_descriptor_list <- c(second_twentyseventh,twentyseventh_GSM,twentyeigth_sixtyseventh)
sample_descriptor_vector <- unlist(sample_descriptor_list)
sample_metadata <- read.table(text = sample_descriptor_vector, sep = " ", colClasses = "character")
colnames(sample_metadata) <- c("sample_no","animal_id","sex","age_days","diet")
sample_metadata$geo_accession <- unlist(sample_id)
sample_title_vector <- unlist(sample_title)
sample_title_df <- read.table(text = sample_title_vector, sep = ",", colClasses = "character")
colnames(sample_title_df) <- c("tissue","array_id","sample_title")
phenoData <- cbind(sample_title_df[,c(1,2)], sample_metadata)
age_days <- as.integer(phenoData$age)
age_weeks <- age_days/7
age_months <- age_weeks/4
phenoData$age_months <- age_months
ordered <- rep(c("6","17","25"),each=19)
ordered <- c(rep("1",times=10),ordered)
couples <- phenoData[,c(3,9)]
sorted_couples <- couples[order(age_months),]
sorted_couples$age_group <- as.integer(ordered)
#
#
#
library(plyr)
#
#
#
sorted_couples_final <- arrange(sorted_couples,sample_no)
phenoData$age_months <- sorted_couples_final$age_group 
males_adlib <- subset(phenoData,sex=="M" & diet=="A", select=c(1:9))
males_adlib_sample_no <- as.integer(males_adlib$sample_no)
females_adlib <- subset(phenoData,sex=="F" & diet=="A", select=c(1:9))
females_adlib_sample_no <- as.integer(females_adlib$sample_no)
adlib_group <- subset(phenoData, diet=="A", select=c(1:9))
adlib_sample_no <- as.integer(adlib_group$sample_no)
males_cr <- subset(phenoData,sex=="M" & diet=="C", select=c(1:9))
males_cr_sample_no <- as.integer(males_cr$sample_no)
females_cr <- subset(phenoData,sex=="F" & diet=="C", select=c(1:9))
females_cr_sample_no <- as.integer(females_cr$sample_no)
cr_group <- subset(phenoData, diet=="A", select=c(1:9))
cr_group_sample_no <- as.integer(cr_group$sample_no)
GPL738_fData <- fData(gsm)
GPL782_fData <- fData(gsm_2)
phenoData_2 <- phenoData
phenoData_2$geo_accession <- array_id_782
phenoData_2$array_id <- rep(" mouse17K-B",times=67)
rownames(phenoData) <- geo_array_id[1:67]
rownames(phenoData_2) <- geo_array_id[68:134]
gse_table_raw <- read.table(gse_matrix, header = FALSE, sep = "\t", blank.lines.skip = TRUE, comment.char = "!", fill = TRUE, skip = 0, nrows = -1)
colnames(gse_table_raw) <- c("ID_REF","VALUE","RAW")
gse_column_array_ids <- rep(geo_array_id,each=8448)
gse_table_raw$array_id <- gse_column_array_ids
gse_A <- gse_table_raw[1:566016,]
gse_B <- gse_table_raw[566017:length(gse_table_raw$VALUE),]
gse_zscore_A <- matrix(gse_A$VALUE,ncol=67,nrow=8448,byrow=FALSE)
gse_zscore_B <- matrix(gse_B$VALUE,ncol=67,nrow=8448,byrow=FALSE)
colnames(gse_zscore_A) <- geo_array_id[1:67]
colnames(gse_zscore_B) <- geo_array_id[68:134]
rownames(gse_zscore_A) <- c(1:8448)
rownames(gse_zscore_B) <- c(8449:16896)
gse_phenoData_mat_A <- new("AnnotatedDataFrame",data=phenoData)
gse_phenoData_mat_B <- new("AnnotatedDataFrame",data=phenoData_2)
gse_featureData_mat_A <- new("AnnotatedDataFrame",data=fData(gsm))
gse_featureData_mat_B <- new("AnnotatedDataFrame",data=fData(gsm_2))
gse_mSet_zscore_A <- ExpressionSet(assayData=gse_zscore_A,phenoData=gse_phenoData_mat_A,featureData=gse_featureData_mat_A,annotation="gpl738")
gse_mSet_zscore_B <- ExpressionSet(assayData=gse_zscore_B,phenoData=gse_phenoData_mat_B,featureData=gse_featureData_mat_B,annotation="gpl782")
gse_matA_preBG_preNorm <- matrix(gse_A$RAW,ncol=67,nrow=8448,byrow=FALSE)
gse_matB_preBG_preNorm <- matrix(gse_B$RAW,ncol=67,nrow=8448,byrow=FALSE)
#
#
#
library(preprocessCore)
#
#
#
gse_matA_preNorm <- rma.background.correct(gse_matA_preBG_preNorm, copy=TRUE)
gse_matB_preNorm <- rma.background.correct(gse_matB_preBG_preNorm, copy=TRUE)
colnames(gse_matA_preBG_preNorm) <- geo_array_id[1:67]
colnames(gse_matB_preBG_preNorm) <- geo_array_id[68:134]
rownames(gse_matA_preBG_preNorm) <- c(1:8448)
rownames(gse_matB_preBG_preNorm) <- c(8449:16896)
colnames(gse_matA_preNorm) <- geo_array_id[1:67]
colnames(gse_matB_preNorm) <- geo_array_id[68:134]
rownames(gse_matA_preNorm) <- c(1:8448)
rownames(gse_matB_preNorm) <- c(8449:16896)
gse_mSet_raw_A <- ExpressionSet(assayData=gse_matA_preBG_preNorm,phenoData=gse_phenoData_mat_A,featureData=gse_featureData_mat_A,annotation="gpl738")
gse_mSet_raw_B <- ExpressionSet(assayData=gse_matB_preBG_preNorm,phenoData=gse_phenoData_mat_B,featureData=gse_featureData_mat_B,annotation="gpl782")
gse_mSet_raw_BG_A <- ExpressionSet(assayData=gse_matA_preNorm,phenoData=gse_phenoData_mat_A,featureData=gse_featureData_mat_A,annotation="gpl738")
gse_mSet_raw_BG_B <- ExpressionSet(assayData=gse_matB_preNorm,phenoData=gse_phenoData_mat_B,featureData=gse_featureData_mat_B,annotation="gpl782")
#
#
#
gse_matA_qnorm <- normalize.quantiles(gse_matA_preNorm, copy=TRUE)
gse_matB_qnorm <- normalize.quantiles(gse_matB_preNorm, copy=TRUE)
gse_matA_qlog2 <- log2(gse_matA_qnorm[,1:67])
gse_matB_qlog2 <- log2(gse_matB_qnorm[,1:67])
colnames(gse_matA_qlog2) <- geo_array_id[1:67]
colnames(gse_matB_qlog2) <- geo_array_id[1:67]
rownames(gse_matA_qlog2) <- c(1:8448)
rownames(gse_matB_qlog2) <- c(8449:16896)
gse_df_qA <- as.data.frame(gse_matA_qlog2)
gse_df_qB <- as.data.frame(gse_matB_qlog2)
colnames(gse_df_qA) <- geo_array_id[1:67]
colnames(gse_df_qB) <- geo_array_id[1:67]
gse_df_qnorm <- rbind(gse_df_qA,gse_df_qB)
colnames(gse_df_qnorm) <- c(1:67)
gse_mSet_qnorm_A <- ExpressionSet(assayData=gse_matA_qlog2,phenoData=gse_phenoData_mat_A,featureData=gse_featureData_mat_A,annotation="gpl738")
gse_mSet_qnorm_B <- ExpressionSet(assayData=gse_matB_qlog2,phenoData=gse_phenoData_mat_B,featureData=gse_featureData_mat_B,annotation="gpl782")
#
#
#
library(arrayQualityMetrics)
library(lattice)
arrayQualityMetrics(expressionset=gse_mSet_zscore_A,outdir="report_z_score_gse_gpl738_age",intgroup=c("age_months", "diet","sex"))
arrayQualityMetrics(expressionset=gse_mSet_zscore_B,outdir="report_z_score_gse_gpl782_age",intgroup=c("age_months", "diet","sex"))
bwplot(VALUE~array_id,data=gse_A)
bwplot(VALUE~array_id,data=gse_A,do.out=FALSE,ylim=c(-3.5,3.5))
#
#
#
arrayQualityMetrics(expressionset=gse_mSet_raw_A,outdir="report_raw_log2_gse_gpl738_age",do.logtransform=TRUE,intgroup=c("age_months", "diet","sex"))
arrayQualityMetrics(expressionset=gse_mSet_raw_B,outdir="report_raw_log2_gse_gpl782_age",do.logtransform=TRUE,intgroup=c("age_months", "diet","sex"))
bwplot(RAW~array_id,data=gse_A)
bwplot(log2(RAW)~array_id,data=gse_A)
bwplot(log2(RAW)~array_id,data=gse_A,do.out=FALSE,ylim=c(3,16))
#
#
#
arrayQualityMetrics(expressionset=gse_mSet_raw_BG_A,outdir="report_rawBG_log2_gse_gpl738_age",do.logtransform=TRUE,intgroup=c("age_months", "diet","sex"))
arrayQualityMetrics(expressionset=gse_mSet_raw_BG_B,outdir="report_rawBG_log2_gse_gpl782_age",do.logtransform=TRUE,intgroup=c("age_months", "diet","sex"))
boxplot(RAW~array_id,data=gse_A)
boxplot(gse_matA_preNorm,use.cols=TRUE,outline=FALSE)
boxplot(RAW~array_id,data=gse_B)
boxplot(gse_matB_preNorm,use.cols=TRUE,outline=FALSE)
#
#
#
arrayQualityMetrics(expressionset=gse_mSet_qnorm_A,outdir="report_qnorm_log2_gse_gpl738_age",intgroup=c("age_months", "diet","sex"))
arrayQualityMetrics(expressionset=gse_mSet_qnorm_B,outdir="report_qnorm_log2_gse_gpl782_age",intgroup=c("age_months", "diet","sex"))
boxplot.matrix(gse_matA_qlog2)
boxplot.matrix(gse_matB_qlog2)