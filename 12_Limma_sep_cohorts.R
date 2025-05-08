library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
library(statmod)

load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")

#Filter on the desired time points 
Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]


# Oslo -------------------------------------------------
Oslo_MoM_exprset <-  Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Cohort %in% "Oslo"]
samp_oslo <- pData(Oslo_MoM_exprset)

#Merge health (preeclamptic or normal) with Time
Timepoint_Group_oslo <- as.factor(paste(Oslo_MoM_exprset$Group, Oslo_MoM_exprset$Timepoint, sep = "."))

#Make a design matrix with timepoint and adjust for age, nulliparity, race and BMI 
design_Oslo <- model.matrix(~0+Timepoint_Group_oslo+as.numeric(Age)+as.factor(Nulliparity)+as.numeric(BMI),samp_oslo)

colnames(design_Oslo) <- c("Control.T2", "Control.T3", "Control.T4", "PE.T2", "PE.T3", "PE.T4", 
                      "Age", "Nulliparity", "BMI")


#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit_oslo <- duplicateCorrelation(Oslo_MoM_exprset, design_Oslo, block = Oslo_MoM_exprset$Cohort_ID)
corfit_oslo$consensus 
#0.541057 Adjusted for Age, Nulliparity, BMI

#Do the limma test, and block for correlations effects from patients
fit_oslo <- lmFit(Oslo_MoM_exprset, design_Oslo, block = Oslo_MoM_exprset$Cohort_ID, correlation = corfit_oslo$consensus)

#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast_oslo <- makeContrasts(T2 = PE.T2 - Control.T2, 
                          T3 = PE.T3 - Control.T3,
                          T4 = PE.T4 - Control.T4, 
                          levels = design_Oslo)

#Compute the moderated t-tests
fit_oslo <- contrasts.fit(fit_oslo, contrast_oslo)
efit_oslo <- eBayes(fit_oslo)

save(fit_oslo, design_Oslo, efit_oslo, file = "../Limma_sep_cohorts/Oslo_fit.RData")


#No cut off values
T2_oslo <- topTable(efit_oslo, coef = "T2", number = "all", adjust.method="BH")
T3_oslo <- topTable(efit_oslo, coef = "T3", number = "all", adjust.method="BH")
T4_oslo <- topTable(efit_oslo, coef = "T4", number = "all", adjust.method="BH")


T2_oslo <- T2_oslo %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3_oslo <- T3_oslo %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T4_oslo <- T4_oslo %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T2_oslo, path = "../Limma_sep_cohorts/T2_oslo.xlsx",
           col_names = TRUE)
write_xlsx(T3_oslo, path = "../Limma_sep_cohorts/T3_oslo.xlsx",
           col_names = TRUE)
write_xlsx(T4_oslo, path = "../Limma_sep_cohorts/T4_oslo.xlsx",
           col_names = TRUE)



#Stanford,  -------------------------------------------------
stanford_MoM_exprset <-  Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Cohort %in% "Stanford"]
samp_stanford <- pData(stanford_MoM_exprset)

#Merge health (preeclamptic or normal) with Time
Timepoint_Group_stanford <- as.factor(paste(stanford_MoM_exprset$Group, stanford_MoM_exprset$Timepoint, sep = "."))

#Make a design matrix with timepoint and adjust for age, nulliparity, race and BMI 
design_stanford <- model.matrix(~0+Timepoint_Group_stanford+as.numeric(Age)+as.numeric(BMI),samp_stanford)


colnames(design_stanford) <- c("Control.T2", "Control.T3", "Control.T4", "PE.T2", "PE.T3", "PE.T4", 
                      "Age", "BMI")



#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit_stanford <- duplicateCorrelation(stanford_MoM_exprset, design_stanford, block = stanford_MoM_exprset$Cohort_ID)
corfit_stanford$consensus 
#0.6770507 Adjusted for Age and BMI

#Do the limma test, and block for correlations effects from patients
fit_stanford <- lmFit(stanford_MoM_exprset, design_stanford, block = stanford_MoM_exprset$Cohort_ID, correlation = corfit_stanford$consensus)

#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast_stanford <- makeContrasts(T2 = PE.T2 - Control.T2, 
                               T3 = PE.T3 - Control.T3,
                               T4 = PE.T4 - Control.T4, 
                               levels = design_stanford)

#Compute the moderated t-tests
fit_stanford <- contrasts.fit(fit_stanford, contrast_stanford)
efit_stanford <- eBayes(fit_stanford)

save(fit_stanford, design_stanford, efit_stanford, file = "../Limma_sep_cohorts/stanford_fit.RData")


#No cut off values
T2_stanford <- topTable(efit_stanford, coef = "T2", number = "all", adjust.method="BH")
T3_stanford <- topTable(efit_stanford, coef = "T3", number = "all", adjust.method="BH")
T4_stanford <- topTable(efit_stanford, coef = "T4", number = "all", adjust.method="BH")


T2_stanford <- T2_stanford %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3_stanford <- T3_stanford %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T4_stanford <- T4_stanford %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T2_stanford, path = "../Limma_sep_cohorts/T2_stanford.xlsx",
           col_names = TRUE)
write_xlsx(T3_stanford, path = "../Limma_sep_cohorts/T3_stanford.xlsx",
           col_names = TRUE)
write_xlsx(T4_stanford, path = "../Limma_sep_cohorts/T4_stanford.xlsx",
           col_names = TRUE)

# detroit -------------------------------------------------
detroit_MoM_exprset <-  Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Cohort %in% "Detroit"]
samp_detroit <- pData(detroit_MoM_exprset)

#Merge health (preeclamptic or normal) with Time
Timepoint_Group_detroit <- as.factor(paste(detroit_MoM_exprset$Group, detroit_MoM_exprset$Timepoint, sep = "."))

#Make a design matrix with timepoint and adjust for age, nulliparity, race and BMI 
design_detroit <- model.matrix(~0+Timepoint_Group_detroit+as.numeric(Age)+as.factor(Nulliparity)+as.numeric(BMI),samp_detroit)


colnames(design_detroit) <- c("Control.T2", "Control.T3", "Control.T4", "PE.T2", "PE.T3", "PE.T4", 
                      "Age", "Nulliparity", "BMI")



#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit_detroit <- duplicateCorrelation(detroit_MoM_exprset, design_detroit, block = detroit_MoM_exprset$Cohort_ID)
corfit_detroit$consensus 
#0.6915134 Adjusted for Age, Nulliparity and BMI

#Do the limma test, and block for correlations effects from patients
fit_detroit <- lmFit(detroit_MoM_exprset, design_detroit, block = detroit_MoM_exprset$Cohort_ID, correlation = corfit_detroit$consensus)

#Make the contrast: Which proteins are deferentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast_detroit <- makeContrasts(T2 = PE.T2 - Control.T2, 
                               T3 = PE.T3 - Control.T3,
                               T4 = PE.T4 - Control.T4, 
                               levels = design_detroit)

#Compute the moderated t-tests
fit_detroit <- contrasts.fit(fit_detroit, contrast_detroit)
efit_detroit <- eBayes(fit_detroit)

save(fit_detroit, design_detroit, efit_detroit, file = "../Limma_sep_cohorts/detroit_fit.RData")


#No cut off values
T2_detroit <- topTable(efit_detroit, coef = "T2", number = "all", adjust.method="BH")
T3_detroit <- topTable(efit_detroit, coef = "T3", number = "all", adjust.method="BH")
T4_detroit <- topTable(efit_detroit, coef = "T4", number = "all", adjust.method="BH")


T2_detroit <- T2_detroit %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3_detroit <- T3_detroit %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T4_detroit <- T4_detroit %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T2_detroit, path = "../Limma_sep_cohorts/T2_detroit.xlsx",
           col_names = TRUE)
write_xlsx(T3_detroit, path = "../Limma_sep_cohorts/T3_detroit.xlsx",
           col_names = TRUE)
write_xlsx(T4_detroit, path = "../Limma_sep_cohorts/T4_detroit.xlsx",
           col_names = TRUE)



