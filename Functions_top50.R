################################ AGE BIN FUNCTIONS ############################################
# Classifies ages by scheme for age binnning 
# @param numeric, single element representing age at death
# @return Character vector
bin_age <- function(x) {
  if (x < 30) {'A'}
  else if (x >= 30 && x < 45) {'B'} 
  else if (x >= 45 && x < 60) {'C'}
  else if (x >= 60 && x < 79) {'D'}
  else {'E'}
}

# Assigns age bin classification
# @param vector, numeric
# @return vector of characters 
set_agebin <- function(x) {
  result <- sapply(x,bin_age)
  result 
}


################################## Extract pData and Bin ######################################
# Extracts phenotype data from GEO Object 
# @param string, GEO object
# @return A list of data frames containing age at death and GEO GSM/Column ID seperated by sex
# @return and a factor column for age binning
# @requires GEOquery, local file containing GSM file
by_tissue <- function(x,f) {
  require(GEOquery)
  tissue <- paste('brain region:', x, sep=' ')
  df <- subset(pData(f), characteristics_ch1.1 == tissue, select = c('characteristics_ch1.5', 'characteristics_ch1.6', 'geo_accession', 'characteristics_ch1.4'))
  df$characteristics_ch1.4 <- gsub('^.{16}', '', df$characteristics_ch1.4)
  df$characteristics_ch1.5 <- gsub('^.{5}', '', df$characteristics_ch1.5)
  df$characteristics_ch1.6 <- as.numeric(gsub('^.{25}', '', df$characteristics_ch1.6))
  df$bin <- as.factor(set_agebin(df$characteristics_ch1.6))
  colnames(df) <- c('sex', 'age', 'geo', 'cause_of_death', 'bin')
  df_male <- df[df$sex == 'Male',]
  df_female <- df[df$sex == 'Female',]
  return(list(df_male, df_female))
}


################################### Calculate T. Drift ########################################
# Calculates transcriptional drift by age bin
# @param data.frame, geo_object
# @return Data frame containing mean drift per sample by age, age bin, cause of death, and 
# @return geo accession value
# @requires GEOquery
transcriptional_drift <- function(tissue, geo_object) {
  require(GEOquery)
  # Extract geo_accession data, males only
  tissue.cc <- row.names(subset(tissue[[1]]))                 
  reference.cc <- row.names(subset(tissue[[1]], bin == 'A'))  
  thirty.cc <- row.names(subset(tissue[[1]], bin == 'B'))     
  fifty.cc <- row.names(subset(tissue[[1]], bin == 'C'))     
  seventy.cc <- row.names(subset(tissue[[1]], bin == 'D'))    
  eighty.cc <- row.names(subset(tissue[[1]], bin == 'E'))    
  
  # Extract expression data, males only
  tissue_mat <- exprs(geo_object)[,tissue.cc]
  
  # Convert to data frame.
  tissue.df <- as.data.frame(tissue_mat)
  
  # Extract top 50% of expressed genes (>= median expression value of young reference)
  tissue.df$mean_expression <- rowMeans(tissue.df[,reference.cc])
  median_expression <- median(tissue.df$mean_expression)
  tissue.df <- subset(tissue.df, mean_expression >= median_expression)
  
  # Calculate expression means for each age bin, no sampling, reference left out of plot
  tissue.df$mean_ref <- rowMeans(tissue.df[,reference.cc])
  tissue.df$mean_thirty <- rowMeans(tissue.df[,thirty.cc])
  tissue.df$mean_fifty <- rowMeans(tissue.df[,fifty.cc])
  tissue.df$mean_seventy <- rowMeans(tissue.df[,seventy.cc])
  tissue.df$mean_eighty <- rowMeans(tissue.df[,eighty.cc])
  
  # Calculate drift values (log fold change == difference), reference left out of plot
  tissue.df$dw_thirty <- tissue.df$mean_thirty - tissue.df$mean_ref
  tissue.df$dw_fifty <- tissue.df$mean_fifty - tissue.df$mean_ref
  tissue.df$dw_seventy <- tissue.df$mean_seventy - tissue.df$mean_ref
  tissue.df$dw_eighty <- tissue.df$mean_eighty - tissue.df$mean_ref
  
  # Return Data frame containing expression values and mean drift by age bin
  tissue.df
}


############################## Calculate Drift Variance #######################################
# Calculates variance for each column in a data frame
# @param data.frame
# @return A numeric vector containing variance of each column in data frame
# @requires plyr
drift_variance <- function(x) {
  require(plyr)
  f <- colwise(var)
  f(x)
}


# Calculates mean transcriptional drift per sample
# @param data.frame, data.frame
# @return Data frame containing mean drift per sample by age, age bin, cause of death, and 
# @return geo accession value
# @requires GEOquery
transcriptional_drift_variance <- function(tissue, df) {
  require(GEOquery)
  # Extract geo_accession data, males only
  tissue.cc <- row.names(subset(tissue[[1]]))                 
  reference.cc <- row.names(subset(tissue[[1]], bin == 'A'))  
  thirty.cc <- row.names(subset(tissue[[1]], bin == 'B'))     
  fifty.cc <- row.names(subset(tissue[[1]], bin == 'C'))     
  seventy.cc <- row.names(subset(tissue[[1]], bin == 'D'))    
  eighty.cc <- row.names(subset(tissue[[1]], bin == 'E')) 
  
  # Extract mean reference expression as a vector
  mean.ref <- df$mean_ref
  
  # Calculate drift for each gene in each sample
  byvar.df <- df[,1:length(tissue.cc)]
  thirty.df <- byvar.df[,thirty.cc] - mean.ref
  fifty.df <- byvar.df[,fifty.cc] - mean.ref
  seventy.df <- byvar.df[,seventy.cc] - mean.ref
  eighty.df <- byvar.df[,eighty.cc] - mean.ref
  
  # Extract total drift variance by sample for each age bin
  thirty.drift <- data.frame(GEO = colnames(thirty.df), VAR = as.numeric(drift_variance(thirty.df)[1,]), AGE = subset(tissue[[1]], bin == 'B')$age, BIN = subset(tissue[[1]], bin == 'B')$bin, COD = subset(tissue[[1]], bin == 'B')$cause_of_death)
    
  fifty.drift <- data.frame(GEO = colnames(fifty.df), VAR = as.numeric(drift_variance(fifty.df)[1,]), AGE = subset(tissue[[1]], bin == 'C')$age, BIN = subset(tissue[[1]], bin == 'C')$bin, COD = subset(tissue[[1]], bin == 'C')$cause_of_death)
    
  seventy.drift <- data.frame(GEO = colnames(seventy.df), VAR = as.numeric(drift_variance(seventy.df)[1,]), AGE = subset(tissue[[1]], bin == 'D')$age, BIN = subset(tissue[[1]], bin == 'D')$bin, COD = subset(tissue[[1]], bin == 'D')$cause_of_death)
    
  eighty.drift <- data.frame(GEO = colnames(eighty.df), VAR = as.numeric(drift_variance(eighty.df)[1,]), AGE = subset(tissue[[1]], bin == 'E')$age, BIN = subset(tissue[[1]], bin == 'E')$bin, COD = subset(tissue[[1]], bin == 'E')$cause_of_death)
    
  drift.bysample <- rbind(thirty.drift, fifty.drift, seventy.drift, eighty.drift)
  
  # Return data frame containing mean drift per sample
  # Factored by age, age bin, cause of death, and geo accession value 
  drift.bysample
}