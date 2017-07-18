# Load Data sets Into R
# Variable descriptions contained in "Creation of Dependent Variable for Q2.docx"
# PC --  C:/Users/Rafael/Documents/DREAM/
# MAC -- ~/DREAM_Prostate_Cancer/
CoreTable_training <- read.csv("~/DREAM_Prostate_Cancer/CoreTable_training.csv", 
                               header=T, na.strings=".", as.is=c("RPT","AGEGRP"), 
                               colClasses=c(STOMACH="factor"))
LabValue_training <- read.csv("~/DREAM_Prostate_Cancer/LabValue_training.csv", 
                              header=T, na.strings=".")
LesionMeasure_training <- read.csv("~/DREAM_Prostate_Cancer/LesionMeasure_training.csv", 
                                   header=T, na.strings=".")
MedHistory_training <- read.csv("~/DREAM_Prostate_Cancer/MedHistory_training.csv", 
                                header=T, na.strings=".")
PriorMed_training <- read.csv("~/DREAM_Prostate_Cancer/PriorMed_training.csv", 
                              header=T, na.strings=".")
VitalSign_training <- read.csv("~/DREAM_Prostate_Cancer/VitalSign_training.csv", 
                               header=T, na.strings=".")
# Data Clean Up for CoreTable_training
#
#
#
# ----- Study Identifiers -----
# RPT == PATIENT ID 
# STUDYID == STUDY SITE IDENTIFIER
#
# ----- Dependent Variables----- 
# LKADT_P == LAST KNOWN ALIVE (Days from Enrollment/Diagnosis?) 
# DEATH == EVENT 
# DISCONT == DISCONTINUED TREATMENT
# ENDTRS_C == REASON TREATMENT ENDED 
# ENDTRT_PC == END TREATMENT DAY (Days from Start of Enrollment/Treatment/Diagnosis?)
#
#
#
# Clinical Variables We Need To Fix
#
# DEATH (''), DISCONT (NA), ENTRT_PC(NA), GLEAS_DX(NA), TSTAG_DX(Missing), HGTBLCAT(NA), 
# WGTBLCAT(NA), SMOKE(Missing,NA), SMOKFREQ(Missing,NA), SMOKSTAT(Missing), 
#
# Clinical Variables Probably Ok To Drop NA/Missing Values
#
# BMI(NA), HEIGHTBL(NA), WEIGHTBL(NA), ECOG_C(NA), TRT2_ID(Missing), TRT3_ID(Missing), 
#
#
# CBC Panel Values With Significant Number of NA Values --- WHY????
#
# LDH(NA, Outliers), TESTO(NA, Outliers), CREACL(NA, Outliers), NA.(NA), MG(NA, Outliers), 
# PHOS(NA, Outliers), ALB(NA), TPRO(NA), RBC(NA), LYM(NA), BUN(NA), CCRC(NA), 
# GLU(NA, Outliers), CREACLCA(NA, Outliers)
#
# CBC Panel Values With < 20 NA Values --- WHY???
#
# ALP(NA, Outliers), ALT(NA, Outliers), AST(NA, Outliers), CA(NA), CREAT(NA, Outliers), 
# HB(NA), NEU(NA, Outliers), PLT(NA, Oultiers), PSA(NA, Outliers), TBILI(NA, Outliers), 
# WBC(NA, Outliers)
#
#
# CBC Panel Values With Low Variance or Small Range (Ignored Outliers Here)
#
# CA, HB, NEU, TESTO, NA., MG, PHOS, RBC, LYM, GLU
#
#
# Need to Fix Factor Levels For Everything Past CREACLCA
#
# c('NO', 'YES') <- c('', 'Y')
# or
#c('NO', 'YES') <- c('', 'YES')
#
#
# 
#
# CONTAINS ALL NA's (1600) or DATA ABSENT -- DROP COLUMNS
# PANCREAS, THYROID, STOMACH
#
#
#
# This is how to fill in missing data/NA
#
# Make copy of original df to manipulate
test_frame <- CoreTable_training
#
table(CoreTable_training$DEATH)
levels(CoreTable_training$DEATH)
# We see "YES" values and empty/blank "" values
# The factor levels are shown as "" and "YES"
# We want to fill in the "" values with default value "NO"
# Need to change levels to "NO" and "YES". 
# The order of appeareance of each factor is important.
# Follow order shown when calling the levels() function.
# This will ensure YES==YES and NO==NO.
levels(test_frame$DEATH) <- c('NO','YES')
#
# Drop Pancreas, Thyroid, Stomach columns
test_frame <- subset(test_frame, select = -c(PANCREAS,THYROID, STOMACH))
#
# Censor Patients with missing discontinue time information.
test_frame <- test_frame[complete.cases(test_frame$ENTRT_PC),]
#
# Censor Patients with miscellaneous end treatment reason 'misce' for ENDTRS_C
test_frame <- subset(test_frame, ENDTRS_C != 'misce')
# Figure out how to drop usued levels, not working right now
levels(droplevels(test_frame$ENDTRS_C))
levels(test_frame$ENDTRS_C) <- c('AE','complete','possible_AE','progression')
#
# Convert DISCONT factor levels to replace 'NA' with 'other' 