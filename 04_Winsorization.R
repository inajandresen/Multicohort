rm(list=ls())
library(tidyverse)
library(Biobase)
library(reshape2)
library(ggforce)
library(pcaMethods)
library(scales)
library(gridExtra)
library(ggpubr)

#Make the function 
winsorize<-function(x)
{
  if(any(x> quantile(x,0.98)*2)){
    x[x>quantile(x,0.98)*2] <- quantile(x,0.98)*2
  }
  return(x)
  
}

#Outlier are defined as any data point above the 2*98th quantile (for that specific protein). 
#The winsorization moves outliers to the 2*98th quantile

# Winsorization Stanford --------------------------------------------------
load("../Data/Merged_cohorts_RFU_exprset.RData")
Stanford_RFU_exprset <- Merged_cohorts_RFU_exprset[, Merged_cohorts_RFU_exprset$Cohort %in% "Stanford"]
Stanford_exprs <- exprs(Stanford_RFU_exprset)

#Are there any out lier?
Stanford_outliers<- apply(Stanford_exprs,1,function(x){any(x> quantile(x,0.98)*2)})
table(Stanford_outliers)
#FALSE  TRUE 
# 940    66 

Stanford_exprs_winsorized <- t(apply(Stanford_exprs,1,winsorize))
par(mfrow= c(2,1))
plot(density(Stanford_exprs["SL000027",]),main="Before")
plot(density(Stanford_exprs_winsorized["SL000027",]),main="After")

#Log2 tranform
Stanford_exprs_log <- log2(Stanford_exprs)
Stanford_exprs_winsorized_log <- log2(Stanford_exprs_winsorized)
par(mfrow= c(2,1))
plot(density(Stanford_exprs_log["SL000027",]),main="Before")
plot(density(Stanford_exprs_winsorized_log["SL000027",]),main="After")


# Winsorisation Detroit ---------------------------------------------------

load("../Data/Merged_cohorts_RFU_exprset.RData")
Detroit_RFU_exprset <- Merged_cohorts_RFU_exprset[, Merged_cohorts_RFU_exprset$Cohort %in% "Detroit"]
Detroit_samp <- pData(Detroit_RFU_exprset)
Detroit_exprs <- exprs(Detroit_RFU_exprset)

#Are there any outliers?
outliers_Detroit <- apply(Detroit_exprs,1,function(x){any(x> quantile(x,0.98)*2)})
table(outliers_Detroit)
#FALSE  TRUE 
# 850   156

#Winsorize the outliers
Detroit_exprs_winsorized <- t(apply(Detroit_exprs,1,winsorize))
outliers_Detroit_win <- apply(Detroit_exprs_winsorized,1,function(x){any(x> quantile(x,0.98)*2)})
table(outliers_Detroit_win)

#Plot one example
par(mfrow= c(2,1))
plot(density(Detroit_exprs["SL000420",]),main="Before")
plot(density(Detroit_exprs_winsorized["SL000420",]),main="After")

#Log2 tranform
Detroit_exprs_log <- log2(Detroit_exprs)
Detroit_exprs_winsorized_log <- log2(Detroit_exprs_winsorized)
par(mfrow= c(2,1))
plot(density(Detroit_exprs_log["SL000420",]),main="Before")
plot(density(Detroit_exprs_winsorized_log["SL000420",]),main="After")

# Winsorisation Oslo ---------------------------------------------------

load("../Data/Merged_cohorts_RFU_exprset.RData")
Oslo_RFU_exprset <- Merged_cohorts_RFU_exprset[, Merged_cohorts_RFU_exprset$Cohort %in% "Oslo"]
Oslo_samp <- pData(Oslo_RFU_exprset)
Oslo_exprs <- exprs(Oslo_RFU_exprset)

#Are there any outliers?
outliers_Oslo <- apply(Oslo_exprs,1,function(x){any(x> quantile(x,0.98)*2)})
table(outliers_Oslo)
#FALSE  TRUE 
# 604   402

#Winsorize the outliers
Oslo_exprs_winsorized <- t(apply(Oslo_exprs,1,winsorize))
outliers_Oslo_win <- apply(Oslo_exprs_winsorized,1,function(x){any(x> quantile(x,0.98)*2)})
table(outliers_Oslo_win)

#Plot one example
par(mfrow= c(2,1))
plot(density(Oslo_exprs["SL016553",]),main="Before")
plot(density(Oslo_exprs_winsorized["SL016553",]),main="After")

#Log2 tranform
Oslo_exprs_log <- log2(Oslo_exprs)
Oslo_exprs_winsorized_log <- log2(Oslo_exprs_winsorized)
par(mfrow= c(2,1))
plot(density(Oslo_exprs_log["SL016553",]),main="Before")
plot(density(Oslo_exprs_winsorized_log["SL016553",]),main="After")


# Merge the Datasets back together ----------------------------------------
samp <- pData(Merged_cohorts_RFU_exprset)
samp <- samp %>% 
  rownames_to_column(var = "rn")
samp <- samp[order(samp$rn),]
samp <- samp %>% 
  column_to_rownames(var = "rn")


prot <- fData(Merged_cohorts_RFU_exprset)

Detroit_exprs_winsorized_log <- Detroit_exprs_winsorized_log %>%
  as.data.frame() %>%
  rownames_to_column(var = "SomaId")

Stanford_exprs_winsorized_log <- Stanford_exprs_winsorized_log %>%
  as.data.frame() %>%
  rownames_to_column(var = "SomaId")

Oslo_exprs_winsorized_log <- Oslo_exprs_winsorized_log %>%
  as.data.frame() %>%
  rownames_to_column(var = "SomaId")

exprs <- merge(Detroit_exprs_winsorized_log, Stanford_exprs_winsorized_log, by = "SomaId")
exprs <- merge(exprs, Oslo_exprs_winsorized_log, by = "SomaId")
exprs <- exprs %>%
  column_to_rownames(var = "SomaId") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "rn")
exprs <- exprs[order(exprs$rn),]
exprs <- exprs %>%
  remove_rownames() %>%
  column_to_rownames(var = "rn") %>%
  t()
  
colnames(exprs)==row.names(samp)

Merged_cohorts_wins_lgRFU_exprset <- ExpressionSet(exprs, phenoData = AnnotatedDataFrame(samp), featureData = AnnotatedDataFrame(prot))
save(Merged_cohorts_wins_lgRFU_exprset, file = "../Data/Merged_cohorts_wins_lgRFU_exprset.RData")
