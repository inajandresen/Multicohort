library(tidyverse)
library(writexl)
library(readxl)
library(Biobase)
library(limma)
#library(ggVennDiagram)
library(ggvenn)
library(statmod)

# Investigate the samples,  -------------------------------------------------
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")

#Filter on the desired time points 
Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]

samp_merged <- pData(Merged_cohorts_wins_Detrended_exprset)

samp_merged$Cohort=as.factor(samp_merged$Cohort)
samp_merged$Nulliparity=as.factor(samp_merged$Nulliparity)
samp_merged$BMI=as.numeric(samp_merged$BMI)
#samp_merged$Race=factor(samp_merged$Race)

prot_merged <- fData(Merged_cohorts_wins_Detrended_exprset)

#How many patients, and how many PE and normal
samp_unique_ids <- samp_merged %>%  #Extract one row per patient
  distinct(Cohort_ID, .keep_all = TRUE)
table(unlist(samp_unique_ids$PEgroups))
#Total: 302 patients 
#PE-late: 124
#Normal: 178

table(unlist(samp_merged$PEgroups))

table(unlist(samp_merged$Timepoint, samp_merged$PEgroups))
#T2: 246
#T3: 335
#T4: 277

#Merge health (preeclamptic or normal) with Time
Timepoint_Group <- as.factor(paste(Merged_cohorts_wins_Detrended_exprset$Group, Merged_cohorts_wins_Detrended_exprset$Timepoint, sep = "."))

#Make a design matrix with timepoint and adjust for age, nulliparity, race and BMI 
design <- model.matrix(~0+Timepoint_Group+Age+Nulliparity+BMI,samp_merged)

colnames(design)<-gsub("Timepoint_Group","",colnames(design)) #will be problems with this later


#See this about blocking: https://mailman.stat.ethz.ch/pipermail/bioconductor/2014-February/057887.html
#Estimate the correlation between time measurements made on the same patient
#Treats patients (ID) as random effect
corfit <- duplicateCorrelation(Merged_cohorts_wins_Detrended_exprset, design, block = Merged_cohorts_wins_Detrended_exprset$Cohort_ID)
corfit$consensus 
#0.654629 Adjusted for Cohort, Age, Nulliparity and BMI 


# Limma t-test ------------------------------------------------------------
#Do the limma test, and block for correlations effects from patients
fit <- lmFit(Merged_cohorts_wins_Detrended_exprset, design, block = Merged_cohorts_wins_Detrended_exprset$Cohort_ID, correlation = corfit$consensus)

#Make the contrast: Which proteins are differentially expressed between preeclamptic and and normal patients 
#at the different time point? 
contrast <- makeContrasts(T2 = PE.T2 - Control.T2, 
                          T3 = PE.T3 - Control.T3,
                          T4 = PE.T4 - Control.T4, 
                          levels = design)

#Compute the moderated t-tests
fit <- contrasts.fit(fit, contrast)
efit <- eBayes(fit)

save(fit, design, efit, file = "../Limma_adj_cohort/fit_adj_cohort.RData")


#Extract the results-------------------------------------------------
load("../Limma_adj_cohort/fit_adj_cohort.RData")

#Cutoff FDR < 0.1 
T2 <- topTable(efit, coef = "T2", number = "all", adjust.method="BH", p.value = 0.1)
T3 <- topTable(efit, coef = "T3", number = "all", adjust.method="BH", p.value = 0.1)
T4 <- topTable(efit, coef = "T4", number = "all", adjust.method="BH", p.value = 0.1)

T2_prot <- T2$EntrezGeneSymbol
T3_prot <- T3$EntrezGeneSymbol
T4_prot <- T4$EntrezGeneSymbol


#Plot Venn diagram for visualization
venn_list <- list("Time interval 1" = T2$SomaId, "Time interval 2" = T3$SomaId, "Time interval 3" = T4$SomaId)

vennplot <- ggvenn(venn_list, 
                   fill_color = c("#A3A500", "#00B0F6", "#E76BF3"),
                   stroke_size = 0.5, 
                   set_name_size = 5,
                   show_percentage = FALSE,
                   text_size = 4.3)

#Total of 62 DE proteins

#No cut off values
T2 <- topTable(efit, coef = "T2", number = "all", adjust.method="BH")
T3 <- topTable(efit, coef = "T3", number = "all", adjust.method="BH")
T4 <- topTable(efit, coef = "T4", number = "all", adjust.method="BH")


T2 <- T2 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T3 <- T3 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

T4 <- T4 %>%
  select(AptName, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)


write_xlsx(T2, path = "../Limma_adj_cohort/T2_adj_cohort.xlsx",
           col_names = TRUE)
write_xlsx(T3, path = "../Limma_adj_cohort/T3_adj_cohort.xlsx",
           col_names = TRUE)
write_xlsx(T4, path = "../Limma_adj_cohort/T4_adj_cohort.xlsx",
           col_names = TRUE)




# How many up- and down regulated -----------------------------------------

T2_all <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T3_all <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T4_all <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")

#Filter on significant proteins, upregulated 
T2_sig_up <- T2_all %>%
  filter(adj.P.Val < 0.1) %>%
  filter(logFC > 0)  
T3_sig_up <- T3_all %>% 
  filter(adj.P.Val < 0.1) %>%
  filter(logFC > 0)  
T4_sig_up <- T4_all %>%
  filter(adj.P.Val < 0.1) %>%
  filter(logFC > 0)  

#Plot Venn diagram for visualization
venn_list_up <- list("Time interval 1" = T2_sig_up$SomaId, "Time interval 2" = T3_sig_up$SomaId, "Time interval 3" = T4_sig_up$SomaId)

vennplot_up <- ggvenn(venn_list_up, 
                   fill_color = c("#A3A500", "#00B0F6", "#E76BF3"),
                   stroke_size = 0.5, 
                   set_name_size = 5,
                   show_percentage = FALSE,
                   text_size = 4.3)

#Total of 40 upregulated

#Filter on significant proteins, downregulated 
T2_sig_down <- T2_all %>%
  filter(adj.P.Val < 0.1) %>%
  filter(logFC < 0)  
T3_sig_down <- T3_all %>% 
  filter(adj.P.Val < 0.1) %>%
  filter(logFC < 0)  
T4_sig_down <- T4_all %>%
  filter(adj.P.Val < 0.1) %>%
  filter(logFC < 0)  

#Plot Venn diagram for visualization
venn_list_down <- list("Time interval 1" = T2_sig_down$SomaId, "Time interval 2" = T3_sig_down$SomaId, "Time interval 3" = T4_sig_down$SomaId)

vennplot_down <- ggvenn(venn_list_down, 
                      fill_color = c("#A3A500", "#00B0F6", "#E76BF3"),
                      stroke_size = 0.5, 
                      set_name_size = 5,
                      show_percentage = FALSE,
                      text_size = 4.3)


