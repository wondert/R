# Load Data sets Into R
# Variable descriptions contained in "Creation of Dependent Variable for Q2.docx"
# PC --  C:/Users/Rafael/Documents/DREAM/
# MAC -- ~/DREAM_Prostate_Cancer/
CoreTable_training <- read.csv("~/DREAM_Prostate_Cancer/CoreTable_training.csv", header=T, na.strings=".", as.is=c("RPT","AGEGRP"), colClasses=c(STOMACH="factor"))
LabValue_training <- read.csv("~/DREAM_Prostate_Cancer/LabValue_training.csv", header=T, na.strings=".")
LesionMeasure_training <- read.csv("~/DREAM_Prostate_Cancer/LesionMeasure_training.csv", header=T, na.strings=".")
MedHistory_training <- read.csv("~/DREAM_Prostate_Cancer/MedHistory_training.csv", header=T, na.strings=".")
PriorMed_training <- read.csv("~/DREAM_Prostate_Cancer/PriorMed_training.csv", header=T, na.strings=".")
VitalSign_training <- read.csv("~/DREAM_Prostate_Cancer/VitalSign_training.csv", header=T, na.strings=".")
#
#
#
# EDA Analysis for CoreTable_training
#
#
#
# ----- Study Identifiers -----
# RPT == PATIENT ID 
# STUDYID == STUDY SITE IDENTIFIER
#
# ----- Dependent Variables----- 
# LKADT_P == LAST KNOWN ALIVE (Days from Enrollment or Diagnosis???) 
# DEATH == EVENT 
# DISCONT == DISCONTINUED TREATMENT
# ENDTRS_C == REASON TREATMENT ENDED 
# ENDTRT_PC == END TREATMENT DAY (Days from Start of Enrollment or Treatment or Diagnosis???)
#
#
# 
# Generate a table containing number of deaths by reason treatment ended to explore relationship
freq_deaths_by_ENDTRSC <- table(deaths=CoreTable_training$DEATH, endtrsc=CoreTable_training$ENDTRS_C)
library(lattice)
barchart(freq_deaths_by_ENDTRSC, stack=FALSE, auto.key=TRUE)
# Patients who died have a greater frequency of "AE", and fare less frequency of "possible_AE"
# Patients who lived are more abundant. Greater frequency of "possible_AE"
# Patients who lived also had greater frequency of "misce" but is likely not relevant
# We should explore relationship between further to include treatment discontinuation
# To accomplish this, I will split groups by DEATH event
pt_death <- subset(CoreTable_training, DEATH=="YES", select=c(1:131))
pt_alive <- subset(CoreTable_training, DEATH=="", select=c(1:131))
freq_DISCONT_by_ENDTRSC_death <- table(deaths=pt_death$DISCONT, endtrsc=pt_death$ENDTRS_C)
freq_DISCONT_by_ENDTRSC_alive <- table(deaths=pt_alive$DISCONT, endtrsc=pt_alive$ENDTRS_C)
barchart(freq_DISCONT_by_ENDTRSC_death, stack=FALSE, auto.key=TRUE)
# Dead patients did not discontinue. No insight there 
# Very high frequency of "AE". This was expected given patient outcome and overall data trends (line 32).
# Lower ~ 100 frequency of "possible_AE" and "progression"
barchart(freq_DISCONT_by_ENDTRSC_alive, stack=FALSE, auto.key=TRUE)
# Significant number of alive patients ~ 100 discontinued treatment. These were "possible_AE" patients.
# Very high frequency of "possible_AE" in alive group. Nearly double the frequency of "AE" alive patient cohort.
# This trend is opposite to what we see in the dead patient pool.
# Frequency of "complete" increased significantly compared to dead patient cohort.
# The "progression" group occurs with similar frequency to dead patient cohort.
#
#
#
# The takeaway from this is that the outcome "possible_AE" is significantly different between the two groups.
#
#
#
# Now we will explore the variable LKADT_P and ENDTRT_PC for total group and alive/dead cohorts
#
#
#