library(GEOquery)
#
#
#
gse_gpl738 <- "GSE7829-GPL738_series_matrix.txt"
gse_gpl782 <- "GSE7829-GPL782_series_matrix.txt"
gsm <- getGEO(filename="GSE7829-GPL738_series_matrix.txt.gz", GSEMatrix=TRUE)
gsm_2 <- getGEO(filename="GSE7829-GPL782_series_matrix.txt.gz", GSEMatrix=TRUE)
gse_matrix <- "GSE7829_family.soft.cleaned_cp.txt"
#
#
#
gse_gpl738_raw <- read.table(gse_gpl738, header = TRUE, sep = "\t", blank.lines.skip = TRUE, comment.char = "!", fill = TRUE, skip = 0, nrows = -1)
gse_gpl782_raw <- read.table(gse_gpl782, header = TRUE, sep = "\t", blank.lines.skip = TRUE, comment.char = "!", fill = TRUE, skip = 0, nrows = -1)
sample_list_numeric <- c(1:67)
array_id_782 <- colnames(gse_gpl782_raw)[2:68]
array_id_738 <- colnames(gse_gpl738_raw)[2:68]
colnames(gse_gpl738_raw) <- c("ID_REF",sample_list_numeric)
colnames(gse_gpl782_raw) <- c("ID_REF",sample_list_numeric)
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