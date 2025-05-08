library(tidyverse)
library(gridExtra)
library(Biobase)
library(ggpubr)
#library(reshape2)
library(reshape)


load("../Data/Merged_cohorts_Detrended_exprset.RData")
Merged_cohorts_Detrended_exprset <- Merged_cohorts_Detrended_exprset[, Merged_cohorts_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]

#Remove dupicated patients 
samp_merged_cohorts <- pData(Merged_cohorts_Detrended_exprset)
samp_merged_cohorts_unique <- samp_merged_cohorts %>%  #Extract one row per patient
  distinct(Cohort_ID, .keep_all = TRUE)

#Create dummy variables
samp_merged_cohorts_unique$Control <- as.factor(samp_merged_cohorts_unique$Control)
samp_merged_cohorts_unique$Cohort_DV <- samp_merged_cohorts_unique$Cohort
samp_merged_cohorts_unique$Cohort_DV <- gsub("Detroit", "0", samp_merged_cohorts_unique$Cohort_DV)
samp_merged_cohorts_unique$Cohort_DV <- gsub("Stanford", "1", samp_merged_cohorts_unique$Cohort_DV)
samp_merged_cohorts_unique$Cohort_DV <- gsub("Oslo", "2", samp_merged_cohorts_unique$Cohort_DV)
samp_merged_cohorts_unique$Cohort_DV <- as.factor(samp_merged_cohorts_unique$Cohort_DV)


#BMI 

#Adjust for cohorts 
BMI_model_adj <- glm(Control ~ BMI+Cohort_DV, family=binomial, data = samp_merged_cohorts_unique)
summary(BMI_model_adj)
#p: 0.000195, p <0.001

#No adjustments
BMI_model <- glm(Control ~ BMI, family=binomial, data = samp_merged_cohorts_unique)
summary(BMI_model)
#p: 5.14e-05, p <0.001



#Age

#Adjust for cohort
Age_model_adj <- glm(Control ~ Age+Cohort_DV, family=binomial, data = samp_merged_cohorts_unique)
summary(Age_model_adj)
#p: 0.341

#No adjustment
Age_model <- glm(Control ~ Age, family=binomial, data = samp_merged_cohorts_unique)
summary(Age_model)
# p = 0.0676


#Nulliparity
samp_merged_cohorts_unique$Nulliparity <- as.factor(samp_merged_cohorts_unique$Nulliparity)

#Adjustment for cohort
Nulli_model_adj <- glm(Control ~ Nulliparity+Cohort_DV, family= "binomial", data = samp_merged_cohorts_unique)
summary(Nulli_model_adj)
# p = 0.0124


#No adjustments
Nulli_model <- glm(Control ~ Nulliparity, family= "binomial", data = samp_merged_cohorts_unique)
summary(Nulli_model)
# p = 0.056603



#Ethnicity
samp_merged_cohorts_unique$Race_DV <- samp_merged_cohorts_unique$Race
samp_merged_cohorts_unique$Race_DV <- gsub("African American", "0", samp_merged_cohorts_unique$Race_DV)
samp_merged_cohorts_unique$Race_DV <- gsub("Caucasian", "1", samp_merged_cohorts_unique$Race_DV)
samp_merged_cohorts_unique$Race_DV <- gsub("Asian", "2", samp_merged_cohorts_unique$Race_DV)
samp_merged_cohorts_unique$Race_DV <- as.factor(samp_merged_cohorts_unique$Race_DV)

#Adjustment for cohort
Race_model_adj <- glm(Control ~ Race_DV+Cohort_DV, family=binomial, data = samp_merged_cohorts_unique)
summary(Race_model_adj)
#Race_DV1, p = 0.457
#Race_DV2, p = 0.981

Race_model_adj0 <- glm(Control ~ Cohort_DV, family=binomial, data = samp_merged_cohorts_unique)
summary(Race_model_adj0)

Race_anov <- anova(Race_model_adj, Race_model_adj0, test = "LRT") # likelihood-ratio test for nested models: test for effect of Race_DV

#How many samples

table(samp_merged_cohorts$Group, samp_merged_cohorts$Timepoint)
