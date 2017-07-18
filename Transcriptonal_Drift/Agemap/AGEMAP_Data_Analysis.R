# Here we will do a sub-analysis for drift with age by sex and normalization method
#
#
#
# The first data set will be derived from Z-scores
#
#
#
# Our first test group will be males.
# Calculate mean z-score by age then determine the z-ratio for each age group as compared to young reference
# Extract z-score data without the ID_REF column of gse7829
male_AL_gse <- gse[,males_adlib_sample_no+1]
males_adlib <- males_adlib[order(males_adlib$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
males_adlib_by_age <- as.integer(males_adlib$sample_no)
m_age_1 <- males_adlib_by_age[1:5]
m_age_6 <- males_adlib_by_age[6:10]
m_age_17 <- males_adlib_by_age[11:15]
m_age_25 <- males_adlib_by_age[16:20]
# Need to make males_adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(male_AL_gse) %in% age_1)
# match(age_1,names(male_AL_gse))
male_AL_gse$m1 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_1])
male_AL_gse$m6 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_6])
male_AL_gse$m17 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_17])
male_AL_gse$m25 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_25])
# Need to calculate z-score difference, SD of z-score difference, and z-ratio between means of z-scores per age group. age_1 mean z-score used as a reference.
male_AL_gse$zdiff_6 <- male_AL_gse$m6 - male_AL_gse$m1
male_AL_gse$zdiff_17 <- male_AL_gse$m17 - male_AL_gse$m1
male_AL_gse$zdiff_25 <- male_AL_gse$m25 - male_AL_gse$m1
male_AL_gse$sd6 <- sd(male_AL_gse$zdiff_6)
male_AL_gse$sd17 <- sd(male_AL_gse$zdiff_17)
male_AL_gse$sd25 <- sd(male_AL_gse$zdiff_25)
male_AL_gse$zratio_6 <- male_AL_gse$zdiff_6/male_AL_gse$sd6
male_AL_gse$zratio_17 <- male_AL_gse$zdiff_17/male_AL_gse$sd17
male_AL_gse$zratio_25 <- male_AL_gse$zdiff_25/male_AL_gse$sd25
# Drift plot of aging thymus visualizing Z-ratio's with age relative to 1 month old animals
boxplot(male_AL_gse$zratio_6,male_AL_gse$zratio_17,male_AL_gse$zratio_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Z Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Male Thymus', outline=FALSE)
#
#
#
# Our second test group will be females.
# Calculate mean z-score by age then determine the z-ratio for each age group as compared to young reference
# Extract z-score data without the ID_REF column of gse7829
female_AL_gse <- gse[,females_adlib_sample_no+1]
females_adlib <- females_adlib[order(females_adlib$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
females_adlib_by_age <- as.integer(females_adlib$sample_no)
f_age_1 <- females_adlib_by_age[1:5]
f_age_6 <- females_adlib_by_age[6:10]
f_age_17 <- females_adlib_by_age[11:15]
f_age_25 <- females_adlib_by_age[16:20]
# Need to make females_adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(female_AL_gse) %in% age_1)
# match(age_1,names(female_AL_gse))
female_AL_gse$m1 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_1])
female_AL_gse$m6 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_6])
female_AL_gse$m17 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_17])
female_AL_gse$m25 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_25])
# Need to calculate z-score difference, SD of z-score difference, and z-ratio between means of z-scores per age group. age_1 mean z-score used as a reference.
female_AL_gse$zdiff_6 <- female_AL_gse$m6 - female_AL_gse$m1
female_AL_gse$zdiff_17 <- female_AL_gse$m17 - female_AL_gse$m1
female_AL_gse$zdiff_25 <- female_AL_gse$m25 - female_AL_gse$m1
female_AL_gse$sd6 <- sd(female_AL_gse$zdiff_6)
female_AL_gse$sd17 <- sd(female_AL_gse$zdiff_17)
female_AL_gse$sd25 <- sd(female_AL_gse$zdiff_25)
female_AL_gse$zratio_6 <- female_AL_gse$zdiff_6/female_AL_gse$sd6
female_AL_gse$zratio_17 <- female_AL_gse$zdiff_17/female_AL_gse$sd17
female_AL_gse$zratio_25 <- female_AL_gse$zdiff_25/female_AL_gse$sd25
# Drift plot of aging thymus visualizing Z-ratio's with age relative to 1 month old animals
boxplot(female_AL_gse$zratio_6,female_AL_gse$zratio_17,female_AL_gse$zratio_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Z Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Female Thymus', outline=FALSE)
#
#
#
# Our third test group will be all ad libitum fed animals.
# Calculate mean z-score by age then determine the z-ratio for each age group as compared to young reference
# Extract z-score data without the ID_REF column of gse7829
all_AL_gse <- gse[,adlib_sample_no+1]
adlib_group <- adlib_group[order(adlib_group$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
adlib_by_age <- as.integer(adlib_group$sample_no)
age_1 <- adlib_by_age[1:5]
age_6 <- adlib_by_age[6:10]
age_17 <- adlib_by_age[11:15]
age_25 <- adlib_by_age[16:20]
# Need to make adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(all_AL_gse) %in% age_1)
# match(age_1,names(all_AL_gse))
all_AL_gse$m1 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_1])
all_AL_gse$m6 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_6])
all_AL_gse$m17 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_17])
all_AL_gse$m25 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_25])
# Need to calculate z-score difference, SD of z-score difference, and z-ratio between means of z-scores per age group. age_1 mean z-score used as a reference.
all_AL_gse$zdiff_6 <- all_AL_gse$m6 - all_AL_gse$m1
all_AL_gse$zdiff_17 <- all_AL_gse$m17 - all_AL_gse$m1
all_AL_gse$zdiff_25 <- all_AL_gse$m25 - all_AL_gse$m1
all_AL_gse$sd6 <- sd(all_AL_gse$zdiff_6)
all_AL_gse$sd17 <- sd(all_AL_gse$zdiff_17)
all_AL_gse$sd25 <- sd(all_AL_gse$zdiff_25)
all_AL_gse$zratio_6 <- all_AL_gse$zdiff_6/all_AL_gse$sd6
all_AL_gse$zratio_17 <- all_AL_gse$zdiff_17/all_AL_gse$sd17
all_AL_gse$zratio_25 <- all_AL_gse$zdiff_25/all_AL_gse$sd25
# Drift plot of aging thymus visualizing Z-ratio's with age relative to 1 month old animals
boxplot(all_AL_gse$zratio_6,all_AL_gse$zratio_17,all_AL_gse$zratio_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Z Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Thymus', outline=FALSE)
#
#
#
# The second data set will be derived from quantile normalized array data
#
#
#
# Here we generate the data as a matrix, then convert to a data frame with identical structure to gse7829
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
# Reset colnames after merge with rbind to characters 1-67
colnames(gse_df_qnorm) <- c(1:67)
# Our first test group will be males.
# Calculate mean expression by age then determine the log fold change for each age group as compared to young reference
# Extract expression data. This data set does not contain the ID_REF column (so no need for the +1).
male_AL_gse <- gse_df_qnorm[,males_adlib_sample_no]
males_adlib <- males_adlib[order(males_adlib$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
males_adlib_by_age <- as.integer(males_adlib$sample_no)
m_age_1 <- males_adlib_by_age[1:5]
m_age_6 <- males_adlib_by_age[6:10]
m_age_17 <- males_adlib_by_age[11:15]
m_age_25 <- males_adlib_by_age[16:20]
# Need to make males_adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(male_AL_gse) %in% age_1)
# match(age_1,names(male_AL_gse))
male_AL_gse$m1 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_1])
male_AL_gse$m6 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_6])
male_AL_gse$m17 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_17])
male_AL_gse$m25 <- rowMeans(male_AL_gse[,names(male_AL_gse) %in% m_age_25])
# Need to calculate log fold change (difference) between means per age group. age_1 mean expression used as a reference.
male_AL_gse$dw_6 <- male_AL_gse$m6 - male_AL_gse$m1
male_AL_gse$dw_17 <- male_AL_gse$m17 - male_AL_gse$m1
male_AL_gse$dw_25 <- male_AL_gse$m25 - male_AL_gse$m1
# Drift plot of aging thymus visualizing log fold change in expression with age relative to 1 month old animals
boxplot(male_AL_gse$dw_6,male_AL_gse$dw_17,male_AL_gse$dw_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Log Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Male Thymus', outline=FALSE)
#
#
#
# Our second test group will be females.
# Calculate mean expression by age then determine the log fold change for each age group as compared to young reference
# Extract expression data. This data set does not contain the ID_REF column (so no need for the +1).
female_AL_gse <- gse_df_qnorm[,females_adlib_sample_no]
females_adlib <- females_adlib[order(females_adlib$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
females_adlib_by_age <- as.integer(females_adlib$sample_no)
f_age_1 <- females_adlib_by_age[1:5]
f_age_6 <- females_adlib_by_age[6:10]
f_age_17 <- females_adlib_by_age[11:15]
f_age_25 <- females_adlib_by_age[16:20]
# Need to make females_adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(female_AL_gse) %in% age_1)
# match(age_1,names(female_AL_gse))
female_AL_gse$m1 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_1])
female_AL_gse$m6 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_6])
female_AL_gse$m17 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_17])
female_AL_gse$m25 <- rowMeans(female_AL_gse[,names(female_AL_gse) %in% f_age_25])
# Need to calculate log fold change (difference) between means per age group. age_1 mean expression used as a reference.
female_AL_gse$dw_6 <- female_AL_gse$m6 - female_AL_gse$m1
female_AL_gse$dw_17 <- female_AL_gse$m17 - female_AL_gse$m1
female_AL_gse$dw_25 <- female_AL_gse$m25 - female_AL_gse$m1
# Drift plot of aging thymus visualizing log fold change in expression with age relative to 1 month old animals
boxplot(female_AL_gse$dw_6,female_AL_gse$dw_17,female_AL_gse$dw_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Log Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Female Thymus', outline=FALSE)
#
#
#
# Our third test group will be all ad libitum fed animals.
# Calculate mean expression by age then determine the log fold change for each age group as compared to young reference
# Extract expression data. This data set does not contain the ID_REF column (so no need for the +1).
all_AL_gse <- gse_df_qnorm[,adlib_sample_no]
adlib_group <- adlib_group[order(adlib_group$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
adlib_by_age <- as.integer(adlib_group$sample_no)
age_1 <- adlib_by_age[1:5]
age_6 <- adlib_by_age[6:10]
age_17 <- adlib_by_age[11:15]
age_25 <- adlib_by_age[16:20]
# Need to make adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
# which(names(all_AL_gse) %in% age_1)
# match(age_1,names(all_AL_gse))
all_AL_gse$m1 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_1])
all_AL_gse$m6 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_6])
all_AL_gse$m17 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_17])
all_AL_gse$m25 <- rowMeans(all_AL_gse[,names(all_AL_gse) %in% age_25])
# Need to calculate log fold change (difference) between means per age group. age_1 mean expression used as a reference.
all_AL_gse$dw_6 <- all_AL_gse$m6 - all_AL_gse$m1
all_AL_gse$dw_17 <- all_AL_gse$m17 - all_AL_gse$m1
all_AL_gse$dw_25 <- all_AL_gse$m25 - all_AL_gse$m1
# Drift plot of aging thymus visualizing Z-ratio's with age relative to 1 month old animals
boxplot(all_AL_gse$dw_6,all_AL_gse$dw_17,all_AL_gse$dw_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Log Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging Thymus', outline=FALSE)
#
#
#
# Now we will try CR animals
#
#
#
# Calculate mean z-score by age then determine the z-ratio for each age group as compared to young reference
# Extract z-score data without the ID_REF column of gse
male_CR_gse <- gse[,males_cr_sample_no+1]
males_cr <- males_cr[order(males_cr$age_months),]
# Use list below, blocked by groups of 5, to get expression value RowMeans of male_AL_gse data frame
males_cr_by_age <- as.integer(males_cr$sample_no)
m_age_1 <- males_cr_by_age[1:5]
m_cr_age_6 <- males_cr_by_age[1:4]
m_cr_age_17 <- males_cr_by_age[5:9]
m_cr_age_25 <- males_cr_by_age[10:12]
# Need to make males_adlib_by_age into integer vector since the values in data frame are characters and the first column of list would be "01", which is not a valid column index
male_CR_gse$m1 <- rowMeans(male_CR_gse[,names(male_CR_gse) %in% m_age_1])
male_CR_gse$m6 <- rowMeans(male_CR_gse[,names(male_CR_gse) %in% m_cr_age_6])
male_CR_gse$m17 <- rowMeans(male_CR_gse[,names(male_CR_gse) %in% m_cr_age_17])
male_CR_gse$m25 <- rowMeans(male_CR_gse[,names(male_CR_gse) %in% m_cr_age_25])
# Check NA values. If any array had an NA value, the $sd column that contain that array will give only NA values
# In this data set there are not enough CR mice. Something is going on here.
# Need to calculate z-score difference, SD of z-score difference, and z-ratio between means of z-scores per age group. age_1 mean z-score used as a reference.
male_CR_gse$zdiff_6 <- male_CR_gse$m6 - male_CR_gse$m1
male_CR_gse$zdiff_17 <- male_CR_gse$m17 - male_CR_gse$m1
male_CR_gse$zdiff_25 <- male_CR_gse$m25 - male_CR_gse$m1
male_CR_gse$sd6 <- sd(male_CR_gse$zdiff_6)
male_CR_gse$sd17 <- sd(male_CR_gse$zdiff_17)
male_CR_gse$sd25 <- sd(male_CR_gse$zdiff_25)
male_CR_gse$zratio_6 <- male_CR_gse$zdiff_6/male_CR_gse$sd6
male_CR_gse$zratio_17 <- male_CR_gse$zdiff_17/male_CR_gse$sd17
male_CR_gse$zratio_25 <- male_CR_gse$zdiff_25/male_CR_gse$sd25
# Drift plot of aging Hippocampus visualizing Z-ratio's with age relative to 1 month old animals
boxplot(male_CR_gse$zratio_6,male_CR_gse$zratio_17,male_CR_gse$zratio_25, names=c('6mo', '17mo', '25mo'), xlab='Mus Musculus', ylab='Z Ratio Expression (Rel. 1mo)', main='Transcriptional Drift of Aging CR Male Cerebral Cortex', outline=FALSE)