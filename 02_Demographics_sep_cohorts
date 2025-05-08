library(tidyverse)
library(gridExtra)
library(Biobase)
library(ggpubr)
library(reshape)

#https://bookdown.org/ndphillips/YaRrr/t-test-t-test.html
#https://www.reneshbedre.com/blog/fisher-exact-test.html?utm_content=cmp-true

load("../Data/Merged_cohorts_Detrended_exprset.RData")
Merged_cohorts_Detrended_exprset <- Merged_cohorts_Detrended_exprset[, Merged_cohorts_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]

samp_merge_cohorts <- pData(Merged_cohorts_Detrended_exprset)
samp_merge_cohorts_unique <- samp_merge_cohorts %>%  #Extract one row per patient
  distinct(Cohort_ID, .keep_all = TRUE)

#By group and cohort
table(samp_merge_cohorts_unique$Cohort,
      samp_merge_cohorts_unique$Group)

#          Control PE
#Detroit       90 76
#Oslo          70 35
#Stanford      18 13

# PE vs Controls Oslo -----------------------------------------------------
load("../Data/Merged_cohorts_Detrended_exprset.RData")
Merged_cohorts_Detrended_exprset <- Merged_cohorts_Detrended_exprset[, Merged_cohorts_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]
Oslo_exprset <- Merged_cohorts_Detrended_exprset[, pData(Merged_cohorts_Detrended_exprset)$Cohort %in% "Oslo"]

samp_Oslo <- pData(Oslo_exprset)
samp_Oslo_unique <- samp_Oslo %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
samp_Oslo_unique_PE <- samp_Oslo_unique[samp_Oslo_unique$Group == "PE", ]
samp_Oslo_unique_Ctrl <- samp_Oslo_unique[samp_Oslo_unique$Group == "Control", ]

#Age
var(samp_Oslo_unique_PE$Age)
var(samp_Oslo_unique_Ctrl$Age)
#Both welch and student t-tests can be used, I guess

t.test(x = samp_Oslo_unique_PE$Age, y = samp_Oslo_unique_Ctrl$Age)
#Welch: p-value = 0.003798 
t.test(x = samp_Oslo_unique_PE$Age, y = samp_Oslo_unique_Ctrl$Age, var.equal = TRUE)
#Student: p-value = 0.002858
mean(samp_Oslo_unique_PE$Age)
sd(samp_Oslo_unique_PE$Age)
mean(samp_Oslo_unique_Ctrl$Age)
sd(samp_Oslo_unique_Ctrl$Age)

#BMI
var(samp_Oslo_unique_PE$BMI)
var(samp_Oslo_unique_Ctrl$BMI)
#Welch should be used

t.test(x = samp_Oslo_unique_PE$BMI, y = samp_Oslo_unique_Ctrl$BMI)
#Welch: p-value = 0.001481
#t.test(x = samp_Oslo_unique_PE$BMI, y = samp_Oslo_unique_Ctrl$BMI, var.equal = TRUE)
#Student: p-value = 0.000163
mean(samp_Oslo_unique_PE$BMI)
sd(samp_Oslo_unique_PE$BMI)
mean(samp_Oslo_unique_Ctrl$BMI)
sd(samp_Oslo_unique_Ctrl$BMI)

#Ethnicity
table(samp_Oslo_unique$Group, samp_Oslo_unique$Race)

#PE, caucasian
(35/35)*100
#Ctrl, caucasian
(70/70)*100


#Nulliparity
table(samp_Oslo_unique$Group, samp_Oslo_unique$Nulliparity)
Nulliparity_table <- table(samp_Oslo_unique$Nulliparity,
                           samp_Oslo_unique$Group)
Nulliparity_table <- as.data.frame(Nulliparity_table)
colnames(Nulliparity_table) <- c("Nulliparity", "Group", "Freq")
Nulliparity_table <- cast(Nulliparity_table, Nulliparity~Group)
Nulliparity_table <- Nulliparity_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Nulliparity")

#Get percentage
#PE, nulliparous
(24/35)*100
#PE, multiparous
(11/35)*100
#Ctrl, nulliparous
(35/70)*100
#Ctrl, multiparous
(35/70)*100

#Can we use a chi squared test? 
chisq.test(Nulliparity_table)$expected
#We can do a chi suqred test 
chi_nulliparity <- chisq.test(Nulliparity_table)
chi_nulliparity

#Fisher's exact test 
fisher_nulliparity <- fisher.test(Nulliparity_table)
fisher_nulliparity

#Number of samples
table(samp_Oslo$Group)
table(samp_Oslo$Group, samp_Oslo$Timepoint)


# PE vs Controls Stanford -------------------------------------------------
load("../Data/Merged_cohorts_Detrended_exprset.RData")
Merged_cohorts_Detrended_exprset <- Merged_cohorts_Detrended_exprset[, Merged_cohorts_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]
Stanford_exprset <- Merged_cohorts_Detrended_exprset[, pData(Merged_cohorts_Detrended_exprset)$Cohort %in% "Stanford"]

samp_stanford <- pData(Stanford_exprset)
samp_stanford_unique <- samp_stanford %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
samp_stanford_unique_PE <- samp_stanford_unique[samp_stanford_unique$Group == "PE", ]
samp_stanford_unique_Ctrl <- samp_stanford_unique[samp_stanford_unique$Group == "Control", ]

#Age
var(samp_stanford_unique_PE$Age)
var(samp_stanford_unique_Ctrl$Age)
#Welch should be used

t.test(x = samp_stanford_unique_PE$Age, y = samp_stanford_unique_Ctrl$Age)
#Welch: p-value = 0.6935
#t.test(x = samp_stanford_unique_PE$Age, y = samp_stanford_unique_Ctrl$Age, var.equal = TRUE)
#Student: p-value = 0.6871
mean(samp_stanford_unique_PE$Age)
sd(samp_stanford_unique_PE$Age)
mean(samp_stanford_unique_Ctrl$Age)
sd(samp_stanford_unique_Ctrl$Age)

#BMI
var(samp_stanford_unique_PE$BMI)
var(samp_stanford_unique_Ctrl$BMI)
#Welch test should be used

t.test(x = samp_stanford_unique_PE$BMI, y = samp_stanford_unique_Ctrl$BMI)
#Welch: p-value = 0.03721
#t.test(x = samp_stanford_unique_PE$BMI, y = samp_stanford_unique_Ctrl$BMI, var.equal = TRUE)
#Student: p-value = 0.01962
mean(samp_stanford_unique_PE$BMI)
sd(samp_stanford_unique_PE$BMI)
mean(samp_stanford_unique_Ctrl$BMI)
sd(samp_stanford_unique_Ctrl$BMI)

#Ethnicity
Race_table <- table(samp_stanford_unique$Group, samp_stanford_unique$Race)
Race_table <- as.data.frame(Race_table)
colnames(Race_table) <- c("Group", "Race", "Freq")
Race_table <- cast(Race_table, Race~Group)
Race_table <- Race_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Race")

#Can we use a chi squared test? 
chisq.test(Race_table)$expected
#We should do a fisher's exact test 

#Do the fisher's exact test 
fisher_race <- fisher.test(Race_table)
fisher_race

#Percentages
#PE, african american
(1/13)*100
#PE, asian
(2/13)*100
#PE, caucasian
(10/13)*100
#Ctrl, caucasian
(18/18)*100



#Nulliparity
Nulliparity_table <- table(samp_stanford_unique$Nulliparity,
                           samp_stanford_unique$Group)

#PE, multiparous
(13/13)*100
#Ctrl, multiparous
(18/18)*100


#Number of samples
table(samp_stanford$Group)
table(samp_stanford$Group, samp_stanford$Timepoint)



# PE vs Controls Detroit --------------------------------------------------
load("../Data/Merged_cohorts_Detrended_exprset.RData")
Merged_cohorts_Detrended_exprset <- Merged_cohorts_Detrended_exprset[, Merged_cohorts_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]
Detroit_exprset <- Merged_cohorts_Detrended_exprset[, pData(Merged_cohorts_Detrended_exprset)$Cohort %in% "Detroit"]

samp_Detroit <- pData(Detroit_exprset)
samp_Detroit_unique <- samp_Detroit %>%  #Extract one row per patient
  distinct(ID, .keep_all = TRUE)
samp_Detroit_unique_PE <- samp_Detroit_unique[samp_Detroit_unique$Group == "PE", ]
samp_Detroit_unique_Ctrl <- samp_Detroit_unique[samp_Detroit_unique$Group == "Control", ]

#Age
var(samp_Detroit_unique_PE$Age)
var(samp_Detroit_unique_Ctrl$Age)
#Welch should be used

t.test(x = samp_Detroit_unique_PE$Age, y = samp_Detroit_unique_Ctrl$Age)
#Welch: p-value = 0.7787
#t.test(x = samp_Detroit_unique_PE$Age, y = samp_Detroit_unique_Ctrl$Age, var.equal = TRUE)
#Student:  p-value = 0.7722

mean(samp_Detroit_unique_PE$Age)
sd(samp_Detroit_unique_PE$Age)
mean(samp_Detroit_unique_Ctrl$Age)
sd(samp_Detroit_unique_Ctrl$Age)

#BMI
var(samp_Detroit_unique_PE$BMI)
var(samp_Detroit_unique_Ctrl$BMI)
#Welch should be used

t.test(x = samp_Detroit_unique_PE$BMI, y = samp_Detroit_unique_Ctrl$BMI)
#Welch: p-value = 0.02447
t.test(x = samp_Detroit_unique_PE$BMI, y = samp_Detroit_unique_Ctrl$BMI, var.equal = TRUE)
#Student: p-value = 0.02215
mean(samp_Detroit_unique_PE$BMI)
sd(samp_Detroit_unique_PE$BMI)
mean(samp_Detroit_unique_Ctrl$BMI)
sd(samp_Detroit_unique_Ctrl$BMI)

#Ethnicity
Race_table <- table(samp_Detroit_unique$Group, samp_Detroit_unique$Race)
Race_table <- as.data.frame(Race_table)
colnames(Race_table) <- c("Group", "Race", "Freq")
Race_table <- cast(Race_table, Race~Group)
Race_table <- Race_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Race")

#Can we use a chi squared test? 
chisq.test(Race_table)$expected
#We should do a fisher's exact test 

#Do the fisher's exact test 
fisher_race <- fisher.test(Race_table)
fisher_race

#Percentages
#PE, african american
(72/76)*100
#PE, caucasian
(4/76)*100
#Ctrl, african american
(84/90)*100
#Ctrl, caucasian
(6/90)*100


#Nulliparity
Nulliparity_table <- table(samp_Detroit_unique$Nulliparity,
                           samp_Detroit_unique$Group)
Nulliparity_table <- as.data.frame(Nulliparity_table)
colnames(Nulliparity_table) <- c("Nulliparity", "Group", "Freq")
Nulliparity_table <- cast(Nulliparity_table, Nulliparity~Group)
Nulliparity_table <- Nulliparity_table %>%
  remove_rownames() %>%
  column_to_rownames(var = "Nulliparity")

#Get percentage
#PE, nulliparous
(32/76)*100
#PE, multiparous
(44/76)*100
#Ctrl, nulliparous
(26/90)*100
#Ctrl, multiparous
(64/90)*100

#Can we use a chi squared test? 
chisq.test(Nulliparity_table)$expected
#We can do a chi suqred test 
chi_nulliparity <- chisq.test(Nulliparity_table)
chi_nulliparity

#Fisher's exact test 
fisher_nulliparity <- fisher.test(Nulliparity_table)
fisher_nulliparity


#Number of samples
table(samp_Detroit$Group)
table(samp_Detroit$Group, samp_Detroit$Timepoint)

