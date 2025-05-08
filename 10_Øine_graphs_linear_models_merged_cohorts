rm(list=ls())
library(splines)
library(mgcv)
library(tidyverse)
library(Biobase)
library(readxl)
library(reshape2)
library(ggforce)



# Significant proteins ----------------------------------------------------
T2 <- read_excel("../Limma_adj_cohort/T2_adj_cohort.xlsx")
T2 <- T2 %>%
  filter(adj.P.Val < 0.1)
T3 <- read_excel("../Limma_adj_cohort/T3_adj_cohort.xlsx")
T3 <- T3 %>%
  filter(adj.P.Val < 0.1)
T4 <- read_excel("../Limma_adj_cohort/T4_adj_cohort.xlsx")
T4 <- T4 %>%
  filter(adj.P.Val < 0.1)

sig_prot <- c(T2$SomaId, T3$SomaId, T4$SomaId)
sig_prot <- unique(sig_prot)

#MMP7, SL000525
#SIGLEC6, SL005217
#PPID, SL007373
#sig_prot <- c("SL000525", "SL005217", "SL007373")


# Filter on significant proteins ------------------------------------------
load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
Merged_cohorts_wins_Detrended_exprset <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% c("T2", "T3", "T4")]
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)

prot <- fData(Merged_cohorts_wins_Detrended_exprset)
prot2 <- prot %>%
  mutate(EntrezGeneSymbol2 = paste0("(", EntrezGeneSymbol,  ")"))
prot2$TargetFullName_EntrezGeneSymbol <- paste(prot2$TargetFullName, prot2$EntrezGeneSymbol2)
prot2 <- prot2 %>% 
  dplyr::select(SomaId, TargetFullName_EntrezGeneSymbol)

swr = function(string, nwrap=35) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)
prot2$TargetFullName_EntrezGeneSymbol <- swr(prot2$TargetFullName_EntrezGeneSymbol)
prot2 <- prot2[sig_prot , ]

Merged_cohorts_exprs <- t(exprs(Merged_cohorts_wins_Detrended_exprset))
Merged_cohorts_exprs <- Merged_cohorts_exprs[, sig_prot]

Merged_cohorts_samp <- pData(Merged_cohorts_wins_Detrended_exprset)
Merged_cohorts_samp <- Merged_cohorts_samp %>%
  dplyr::select(ID, Group, GAWeeks, Cohort)

Merged_cohorts_exprs_samp <- merge(Merged_cohorts_samp, Merged_cohorts_exprs, by = "row.names")
Merged_cohorts_exprs_samp <- Merged_cohorts_exprs_samp %>%
  column_to_rownames(var = "Row.names")
Merged_cohorts_exprs_samp$ID <- as.factor(Merged_cohorts_exprs_samp$ID)
Merged_cohorts_exprs_samp$Group <- as.factor(Merged_cohorts_exprs_samp$Group)

#Control samples
Control_df <- Merged_cohorts_exprs_samp %>%
  filter(Group == "Control") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

#LOPE samples
LOPE_df <- Merged_cohorts_exprs_samp %>%
  filter(Group == "LOPE") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)


# Fit a model for the Controls --------------------------------------

Pred_df_con <- NULL
protnames <- colnames(Control_df[6:ncol(Control_df)])

for (prot in protnames) {
  
  Control_df$Y=Control_df[,prot]
  fit3_control <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                      data=Control_df,
                      method="ML")
  
  fit4_control <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                      data=Control_df,
                      method="ML")
  
  fit5_control <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                      data=Control_df,
                      method="ML")
  
  fit6_control <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                      data=Control_df,
                      method="ML")
  
  fit_control=list(fit3_control,fit4_control,fit5_control,fit6_control)[[which.max(c(summary(fit3_control)$r.sq,
                                                                                     summary(fit4_control)$r.sq,
                                                                                     summary(fit5_control)$r.sq, 
                                                                                     summary(fit6_control)$r.sq))]]
  
  Control_pred <- as.data.frame(Control_df$GAWeeks)
  colnames(Control_pred) <- "GAWeeks"
  Control_pred$pred <- predict.gam(fit_control,Control_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Control_pred$SomaId=prot
  Pred_df_con=rbind(Pred_df_con,Control_pred)
}

Pred_df_con <- merge(Pred_df_con, prot2, by = "SomaId")
Pred_df_con$Group <- "Control"
save(Control_df, Pred_df_con, file = "../Line_graphs_DE_proteins/Pred_df_controls_DE_proteins.RData")

# Fit a model for the Preeclamptic --------------------------------------

Pred_df_LOPE <- NULL
protnames <- colnames(LOPE_df[6:ncol(LOPE_df)])

for (prot in protnames) {
  
  LOPE_df$Y=LOPE_df[,prot]
  fit3_LOPE <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                   data=LOPE_df,
                   method="ML")
  
  fit4_LOPE <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                   data=LOPE_df,
                   method="ML")
  
  fit5_LOPE <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                   data=LOPE_df,
                   method="ML")
  
  fit6_LOPE <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                   data=LOPE_df,
                   method="ML")
  
  fit_LOPE=list(fit3_LOPE,fit4_LOPE,fit5_LOPE,fit6_LOPE)[[which.max(c(summary(fit3_LOPE)$r.sq,
                                                                      summary(fit4_LOPE)$r.sq,
                                                                      summary(fit5_LOPE)$r.sq, 
                                                                      summary(fit6_LOPE)$r.sq))]]
  
  LOPE_pred <- as.data.frame(LOPE_df$GAWeeks)
  colnames(LOPE_pred) <- "GAWeeks"
  LOPE_pred$pred <- predict.gam(fit_LOPE,LOPE_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  LOPE_pred$SomaId=prot
  Pred_df_LOPE=rbind(Pred_df_LOPE,LOPE_pred)
}

Pred_df_LOPE <- merge(Pred_df_LOPE, prot2, by = "SomaId")
Pred_df_LOPE$Group <- "LOPE"
save(LOPE_df, Pred_df_LOPE, file = "../Line_graphs_DE_proteins/Pred_df_LOPE_DE_proteins.RData")


# Plot all DE proteins --------------------------------------------------------------------
load("../Line_graphs_DE_proteins/Pred_df_controls_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_df_LOPE_DE_proteins.RData")


All_df <- rbind(LOPE_df, Control_df)

All_df_long <- All_df %>%
  dplyr::select(-ID, -Group, -GAWeeks, -Cohort, -Y)
All_df_long <- melt(All_df_long, id = "SampleID")
colnames(All_df_long) <- c("SampleID", "SomaId", "MoM")

All_df_long <- merge(All_df_long, prot2, by = "SomaId")
All_df_samp <- All_df %>%
  dplyr::select(SampleID, ID, Group, GAWeeks, Cohort)
All_df_long <- merge(All_df_samp, All_df_long, by = "SampleID")


#Separating cohorts
pdf("../Line_graphs_DE_proteins/Lineplots_DE_proteins.pdf",
    width = 10, height = 10)

for(i in 1:9) {
  print(ggplot(data = All_df_long, aes(x = GAWeeks, y = MoM)) +
          geom_point(aes(shape = Group, color = Group), alpha = 0.6) +
          scale_color_manual(values = c("grey74", "grey40")) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
          geom_line(data = Pred_df_con, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "Black", size=1) +
          geom_line(data = Pred_df_LOPE, aes(x = GAWeeks, y = pred), linetype = "solid", color = "Black", size=1) + 
          theme_classic() +
          theme(legend.title = element_blank(),
                legend.position = "right",
                strip.text = element_text(size = 13),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 12)) +
          xlab("Gestational weeks") +
          #facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 2, nrow = 3, scales = "free", page = i))
          facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 2, nrow = 3, page = i))
}
dev.off()






