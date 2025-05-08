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


# Fit a model for the Detroit Controls --------------------------------------

#Control samples
Detroit_control_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Detroit") %>%
  filter(Group == "Control") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Detroit_control_df <- NULL
protnames <- colnames(Detroit_control_df[6:ncol(Detroit_control_df)])

for (prot in protnames) {

Detroit_control_df$Y=Detroit_control_df[,prot]

fit3_Detroit_control <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                            data=Detroit_control_df,
                            method="ML")

fit4_Detroit_control <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                            data=Detroit_control_df,
                            method="ML")

fit5_Detroit_control <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                            data=Detroit_control_df,
                            method="ML")

fit6_Detroit_control <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                            data=Detroit_control_df,
                            method="ML")

fit_Detroit_control=list(fit3_Detroit_control,fit4_Detroit_control,
                         fit5_Detroit_control,fit6_Detroit_control)[[which.max(c(summary(fit3_Detroit_control)$r.sq,
                                                                                 summary(fit4_Detroit_control)$r.sq,
                                                                                 summary(fit5_Detroit_control)$r.sq, 
                                                                                 summary(fit6_Detroit_control)$r.sq))]]

Detroit_control_df_pred <- as.data.frame(Detroit_control_df$GAWeeks)
colnames(Detroit_control_df_pred) <- "GAWeeks"
Detroit_control_df_pred$pred <- predict.gam(fit_Detroit_control,Detroit_control_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
Detroit_control_df_pred$SomaId=prot
Pred_Detroit_control_df=rbind(Pred_Detroit_control_df,Detroit_control_df_pred)

}

Pred_Detroit_control_df <- merge(Pred_Detroit_control_df, prot2, by = "SomaId")
save(Detroit_control_df, Pred_Detroit_control_df, file = "../Line_graphs_DE_proteins/Pred_Detroit_control_DE_proteins.RData")

# Fit a model for the Detroit LOPEs --------------------------------------

#LOPE samples
Detroit_LOPE_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Detroit") %>%
  filter(Group == "LOPE") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Detroit_LOPE_df <- NULL
protnames <- colnames(Detroit_LOPE_df[6:ncol(Detroit_LOPE_df)])

for (prot in protnames) {
  
  Detroit_LOPE_df$Y=Detroit_LOPE_df[,prot]
  
  fit3_Detroit_LOPE <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                              data=Detroit_LOPE_df,
                              method="ML")
  
  fit4_Detroit_LOPE <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                              data=Detroit_LOPE_df,
                              method="ML")
  
  fit5_Detroit_LOPE <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                              data=Detroit_LOPE_df,
                              method="ML")
  
  fit6_Detroit_LOPE <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                              data=Detroit_LOPE_df,
                              method="ML")
  
  fit_Detroit_LOPE=list(fit3_Detroit_LOPE,fit4_Detroit_LOPE,
                           fit5_Detroit_LOPE,fit6_Detroit_LOPE)[[which.max(c(summary(fit3_Detroit_LOPE)$r.sq,
                                                                                   summary(fit4_Detroit_LOPE)$r.sq,
                                                                                   summary(fit5_Detroit_LOPE)$r.sq, 
                                                                                   summary(fit6_Detroit_LOPE)$r.sq))]]
  
  Detroit_LOPE_df_pred <- as.data.frame(Detroit_LOPE_df$GAWeeks)
  colnames(Detroit_LOPE_df_pred) <- "GAWeeks"
  Detroit_LOPE_df_pred$pred <- predict.gam(fit_Detroit_LOPE,Detroit_LOPE_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Detroit_LOPE_df_pred$SomaId=prot
  Pred_Detroit_LOPE_df=rbind(Pred_Detroit_LOPE_df,Detroit_LOPE_df_pred)
   
}

Pred_Detroit_LOPE_df <- merge(Pred_Detroit_LOPE_df, prot2, by = "SomaId")
save(Detroit_LOPE_df, Pred_Detroit_LOPE_df, file = "../Line_graphs_DE_proteins/Pred_Detroit_LOPE_DE_proteins.RData")


# Fit a model for the Oslo Controls --------------------------------------

#Control samples
Oslo_control_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Oslo") %>%
  filter(Group == "Control") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Oslo_control_df <- NULL
protnames <- colnames(Oslo_control_df[6:ncol(Oslo_control_df)])

for (prot in protnames) {
  
  Oslo_control_df$Y=Oslo_control_df[,prot]
  
  fit3_Oslo_control <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                              data=Oslo_control_df,
                              method="ML")
  
  fit4_Oslo_control <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                              data=Oslo_control_df,
                              method="ML")
  
  fit5_Oslo_control <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                              data=Oslo_control_df,
                              method="ML")
  
  fit6_Oslo_control <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                              data=Oslo_control_df,
                              method="ML")
  
  fit_Oslo_control=list(fit3_Oslo_control,fit4_Oslo_control,
                           fit5_Oslo_control,fit6_Oslo_control)[[which.max(c(summary(fit3_Oslo_control)$r.sq,
                                                                                   summary(fit4_Oslo_control)$r.sq,
                                                                                   summary(fit5_Oslo_control)$r.sq, 
                                                                                   summary(fit6_Oslo_control)$r.sq))]]
  
  Oslo_control_df_pred <- as.data.frame(Oslo_control_df$GAWeeks)
  colnames(Oslo_control_df_pred) <- "GAWeeks"
  Oslo_control_df_pred$pred <- predict.gam(fit_Oslo_control,Oslo_control_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Oslo_control_df_pred$SomaId=prot
  Pred_Oslo_control_df=rbind(Pred_Oslo_control_df,Oslo_control_df_pred)
  
}

Pred_Oslo_control_df <- merge(Pred_Oslo_control_df, prot2, by = "SomaId")
save(Oslo_control_df, Pred_Oslo_control_df, file = "../Line_graphs_DE_proteins/Pred_Oslo_control_DE_proteins.RData")

# Fit a model for the Oslo LOPEs --------------------------------------

#LOPE samples
Oslo_LOPE_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Oslo") %>%
  filter(Group == "LOPE") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Oslo_LOPE_df <- NULL
protnames <- colnames(Oslo_LOPE_df[6:ncol(Oslo_LOPE_df)])

for (prot in protnames) {
  
  Oslo_LOPE_df$Y=Oslo_LOPE_df[,prot]
  
  fit3_Oslo_LOPE <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                           data=Oslo_LOPE_df,
                           method="ML")
  
  fit4_Oslo_LOPE <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                           data=Oslo_LOPE_df,
                           method="ML")
  
  fit5_Oslo_LOPE <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                           data=Oslo_LOPE_df,
                           method="ML")
  
  fit6_Oslo_LOPE <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                           data=Oslo_LOPE_df,
                           method="ML")
  
  fit_Oslo_LOPE=list(fit3_Oslo_LOPE,fit4_Oslo_LOPE,
                        fit5_Oslo_LOPE,fit6_Oslo_LOPE)[[which.max(c(summary(fit3_Oslo_LOPE)$r.sq,
                                                                          summary(fit4_Oslo_LOPE)$r.sq,
                                                                          summary(fit5_Oslo_LOPE)$r.sq, 
                                                                          summary(fit6_Oslo_LOPE)$r.sq))]]
  
  Oslo_LOPE_df_pred <- as.data.frame(Oslo_LOPE_df$GAWeeks)
  colnames(Oslo_LOPE_df_pred) <- "GAWeeks"
  Oslo_LOPE_df_pred$pred <- predict.gam(fit_Oslo_LOPE,Oslo_LOPE_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Oslo_LOPE_df_pred$SomaId=prot
  Pred_Oslo_LOPE_df=rbind(Pred_Oslo_LOPE_df,Oslo_LOPE_df_pred)
  
}

Pred_Oslo_LOPE_df <- merge(Pred_Oslo_LOPE_df, prot2, by = "SomaId")
save(Oslo_LOPE_df, Pred_Oslo_LOPE_df, file = "../Line_graphs_DE_proteins/Pred_Oslo_LOPE_DE_proteins.RData")

# Fit a model for the Stanford Controls --------------------------------------

#Control samples
Stanford_control_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Stanford") %>%
  filter(Group == "Control") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Stanford_control_df <- NULL
protnames <- colnames(Stanford_control_df[6:ncol(Stanford_control_df)])

for (prot in protnames) {
  
  Stanford_control_df$Y=Stanford_control_df[,prot]
  
  fit3_Stanford_control <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                              data=Stanford_control_df,
                              method="ML")
  
  fit4_Stanford_control <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                              data=Stanford_control_df,
                              method="ML")
  
  fit5_Stanford_control <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                              data=Stanford_control_df,
                              method="ML")
  
  fit6_Stanford_control <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                              data=Stanford_control_df,
                              method="ML")
  
  fit_Stanford_control=list(fit3_Stanford_control,fit4_Stanford_control,
                           fit5_Stanford_control,fit6_Stanford_control)[[which.max(c(summary(fit3_Stanford_control)$r.sq,
                                                                                   summary(fit4_Stanford_control)$r.sq,
                                                                                   summary(fit5_Stanford_control)$r.sq, 
                                                                                   summary(fit6_Stanford_control)$r.sq))]]
  
  Stanford_control_df_pred <- as.data.frame(Stanford_control_df$GAWeeks)
  colnames(Stanford_control_df_pred) <- "GAWeeks"
  Stanford_control_df_pred$pred <- predict.gam(fit_Stanford_control,Stanford_control_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Stanford_control_df_pred$SomaId=prot
  Pred_Stanford_control_df=rbind(Pred_Stanford_control_df,Stanford_control_df_pred)
  
}

Pred_Stanford_control_df <- merge(Pred_Stanford_control_df, prot2, by = "SomaId")
save(Stanford_control_df, Pred_Stanford_control_df, file = "../Line_graphs_DE_proteins/Pred_Stanford_control_DE_proteins.RData")


# Fit a model for the Stanford LOPEs --------------------------------------

#LOPE samples
Stanford_LOPE_df <- Merged_cohorts_exprs_samp %>%
  filter(Cohort == "Stanford") %>%
  filter(Group == "LOPE") %>%
  rownames_to_column(var = "SampleID") %>%
  arrange(GAWeeks)

Pred_Stanford_LOPE_df <- NULL
protnames <- colnames(Stanford_LOPE_df[6:ncol(Stanford_LOPE_df)])

for (prot in protnames) {
  
  Stanford_LOPE_df$Y=Stanford_LOPE_df[,prot]
  
  fit3_Stanford_LOPE <- gam(Y~s(GAWeeks,k=3)+s(ID,bs = 're'),
                           data=Stanford_LOPE_df,
                           method="ML")
  
  fit4_Stanford_LOPE <- gam(Y~s(GAWeeks,k=4)+s(ID,bs = 're'),
                           data=Stanford_LOPE_df,
                           method="ML")
  
  fit5_Stanford_LOPE <- gam(Y~s(GAWeeks,k=5)+s(ID,bs = 're'),
                           data=Stanford_LOPE_df,
                           method="ML")
  
  fit6_Stanford_LOPE <- gam(Y~s(GAWeeks,k=6)+s(ID,bs = 're'),
                           data=Stanford_LOPE_df,
                           method="ML")
  
  fit_Stanford_LOPE=list(fit3_Stanford_LOPE,fit4_Stanford_LOPE,
                        fit5_Stanford_LOPE,fit6_Stanford_LOPE)[[which.max(c(summary(fit3_Stanford_LOPE)$r.sq,
                                                                          summary(fit4_Stanford_LOPE)$r.sq,
                                                                          summary(fit5_Stanford_LOPE)$r.sq, 
                                                                          summary(fit6_Stanford_LOPE)$r.sq))]]
  
  Stanford_LOPE_df_pred <- as.data.frame(Stanford_LOPE_df$GAWeeks)
  colnames(Stanford_LOPE_df_pred) <- "GAWeeks"
  Stanford_LOPE_df_pred$pred <- predict.gam(fit_Stanford_LOPE,Stanford_LOPE_df_pred,exclude="s(ID)",newdata.guaranteed=TRUE)
  Stanford_LOPE_df_pred$SomaId=prot
  Pred_Stanford_LOPE_df=rbind(Pred_Stanford_LOPE_df,Stanford_LOPE_df_pred)
 
}

Pred_Stanford_LOPE_df <- merge(Pred_Stanford_LOPE_df, prot2, by = "SomaId")
save(Stanford_LOPE_df, Pred_Stanford_LOPE_df, file = "../Line_graphs_DE_proteins/Pred_Stanford_LOPE_DE_proteins.RData")


# Merge data and plot -----------------------------------------------------
load("../Line_graphs_DE_proteins/Pred_Stanford_LOPE_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_Stanford_control_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_Oslo_LOPE_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_Oslo_control_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_Detroit_LOPE_DE_proteins.RData")
load("../Line_graphs_DE_proteins/Pred_Detroit_control_DE_proteins.RData")


All_df <- rbind(Detroit_control_df, Detroit_LOPE_df)
All_df <- rbind(All_df, Oslo_control_df)
All_df <- rbind(All_df, Oslo_LOPE_df)
All_df <- rbind(All_df, Stanford_control_df)
All_df <- rbind(All_df, Stanford_LOPE_df)

All_df_long <- All_df %>%
  #dplyr::select(-ID, -Group, -GAWeeks, -Cohort, -Y)
  dplyr::select(-ID, -Group, -GAWeeks, -Cohort)
All_df_long <- melt(All_df_long, id = "SampleID")
colnames(All_df_long) <- c("SampleID", "SomaId", "MoM")

All_df_long <- merge(All_df_long, prot2, by = "SomaId")
All_df_samp <- All_df %>%
  dplyr::select(SampleID, ID, Group, GAWeeks, Cohort)
All_df_long <- merge(All_df_samp, All_df_long, by = "SampleID")


plot <- ggplot(data = All_df_long, aes(x = GAWeeks, y = MoM)) +
  geom_point(aes(shape = Group, color = Group), alpha = 0.6) +
  scale_color_manual(values = c("grey74", "grey40")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_line(data = Pred_Detroit_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#F8766D", size=1) +
  geom_line(data = Pred_Detroit_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#F8766D", size=1) + 
  geom_line(data = Pred_Oslo_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#00BA38", size=1) +
  geom_line(data = Pred_Oslo_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#00BA38", size=1) + 
  geom_line(data = Pred_Stanford_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#619CFF", size=1) +
  geom_line(data = Pred_Stanford_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#619CFF", size=1) + 
  theme_classic() +
  ylab("logMoM") +
  xlab("Weeks") +
  theme(legend.title = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12)) +
  #facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 2, nrow = 3, scales = "free", page = i))
  #facet_wrap_paginate(~TargetFullName_EntrezGeneSymbol, ncol = 4, nrow = 7, page = 2)
  facet_wrap(~TargetFullName_EntrezGeneSymbol, ncol = 5, scales = "free")

ggsave(filename= "Line_plots_DE_proteins_sep_cohorts.png",
       plot = plot,
       device = "png",
       path = "../Line_graphs_DE_proteins/",
       width = 60,
       height = 90,
       units = "cm")




#Separating cohorts
pdf("../Line_graphs_DE_proteins/Lineplots_DE_proteins_sep_cohorts.pdf",
    width = 10, height = 10)

for(i in 1:9) {
  print(ggplot(data = All_df_long, aes(x = GAWeeks, y = MoM)) +
          geom_point(aes(shape = Group, color = Group), alpha = 0.6) +
          scale_color_manual(values = c("grey74", "grey40")) +
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
          geom_line(data = Pred_Detroit_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#F8766D", size=1) +
          geom_line(data = Pred_Detroit_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#F8766D", size=1) + 
          geom_line(data = Pred_Oslo_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#00BA38", size=1) +
          geom_line(data = Pred_Oslo_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#00BA38", size=1) + 
          geom_line(data = Pred_Stanford_control_df, aes(x = GAWeeks, y = pred), linetype = "dashed", color = "#619CFF", size=1) +
          geom_line(data = Pred_Stanford_LOPE_df, aes(x = GAWeeks, y = pred), linetype = "solid", color = "#619CFF", size=1) + 
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



