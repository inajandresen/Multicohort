rm(list=ls())
library(pROC)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpattern)
library(ggpubr)

load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)

Merged_cohorts_wins_Detrended_exprset_T2 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T2"]
samp_T2 <- pData(Merged_cohorts_wins_Detrended_exprset_T2)
samp_T2 <- samp_T2 %>% dplyr::select(Cohort, Group)
Merged_cohorts_wins_Detrended_exprset_T3 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T3"]
samp_T3 <- pData(Merged_cohorts_wins_Detrended_exprset_T3)
samp_T3 <- samp_T3 %>% dplyr::select(Cohort, Group)
Merged_cohorts_wins_Detrended_exprset_T4 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T4"]
samp_T4 <- pData(Merged_cohorts_wins_Detrended_exprset_T4)
samp_T4 <- samp_T4 %>% dplyr::select(Cohort, Group)


# Elastic net Time interval 1 -------------------------------------------------------------

#Time interval 1 (time point 2)
load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T2_EN.RData")
T2_EN_freq <- freq
T2_EN_pile <- pile

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}
T2_EN_a2 <- T2_EN_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")

T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_RC_dt=roc(response=T2_EN_a[grep("Detroit", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Detroit", rownames(T2_EN_a)),2],direction="<")
T2_EN_dt_auc <- round(T2_EN_RC_dt$auc, digits = 2)
T2_EN_dt_ci <- ci(T2_EN_RC_dt, of = "auc")
T2_EN_dt_ci_txt <- "95% CI: 0.54-0.75"

T2_EN_a_dt <- T2_EN_a2[T2_EN_a2$Cohort %in% "Detroit",]

T2_EN_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_a_dt[,2]>=pr,1,0))), factor((T2_EN_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_a_dt[,2]>pr,1,0))), factor(as.character(T2_EN_a_dt[,1])),positive="1")
  T2_EN_tabout_dt=rbind(T2_EN_tabout_dt,c(pr,aa$byClass[1:4]))
}

T2_EN_RC_st=roc(response=T2_EN_a[grep("Stanford", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Stanford", rownames(T2_EN_a)),2],direction="<")
T2_EN_st_auc <- round(T2_EN_RC_st$auc, digits = 2)
T2_EN_st_ci <- ci(T2_EN_RC_st, of = "auc")
T2_EN_st_ci_txt <- "95% CI: 0.21-0.66"

T2_EN_a_st <- T2_EN_a2[T2_EN_a2$Cohort %in% "Stanford",]

T2_EN_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_a_st[,2]>=pr,1,0))), factor((T2_EN_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_a_st[,2]>pr,1,0))), factor(as.character(T2_EN_a_st[,1])),positive="1")
  T2_EN_tabout_st=rbind(T2_EN_tabout_st,c(pr,aa$byClass[1:4]))
}


T2_EN_RC_osl=roc(response=T2_EN_a[grep("Oslo", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Oslo", rownames(T2_EN_a)),2],direction="<")
T2_EN_osl_auc <- round(T2_EN_RC_osl$auc, digits = 2)
T2_EN_osl_ci <- ci(T2_EN_RC_osl, of = "auc")
T2_EN_osl_ci_txt <- "95% CI: 0.42-0.67"

T2_EN_a_osl <- T2_EN_a2[T2_EN_a2$Cohort %in% "Oslo",]

T2_EN_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_a_osl[,2]>=pr,1,0))), factor((T2_EN_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_a_osl[,2]>pr,1,0))), factor(as.character(T2_EN_a_osl[,1])),positive="1")
  T2_EN_tabout_osl=rbind(T2_EN_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
EN_ROC_plot_T2 <- ggroc(list("Detroit" = T2_EN_RC_dt,
                             "Stanford" = T2_EN_RC_st,
                             "Oslo" = T2_EN_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                                           color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_EN_dt_auc, ', ', T2_EN_dt_ci_txt), 
                    paste0('AUC = ', T2_EN_st_auc, ', ', T2_EN_st_ci_txt), 
                    paste0('AUC = ', T2_EN_osl_auc, ', ', T2_EN_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T2_EN_a_mod <- as.data.frame(T2_EN_a)
T2_EN_a_mod <- merge(T2_EN_a_mod, samp_T2, by = "row.names")

T2_EN_box <- ggplot(T2_EN_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



# Elastic net Time interval 2-------------------------------------------------------------

load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T3_EN.RData")
T3_EN_freq <- freq
T3_EN_pile <- pile

T3_EN_a=NULL
for (i in 1:length(T3_EN_pile)){
  T3_EN_a=rbind(T3_EN_a,T3_EN_pile[[i]]$outmat)
}
T3_EN_a2 <- T3_EN_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")


T3_EN_RC=roc(response=T3_EN_a[,1],predictor=T3_EN_a[,2],direction="<")
T3_EN_RC_dt=roc(response=T3_EN_a[grep("Detroit", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Detroit", rownames(T3_EN_a)),2],direction="<")
T3_EN_dt_auc <- round(T3_EN_RC_dt$auc, digits = 2)
T3_EN_dt_ci <- ci(T3_EN_RC_dt, of = "auc")
T3_EN_dt_ci_txt <- "95% CI: 0.49-0.64"

T3_EN_a_dt <- T3_EN_a2[T3_EN_a2$Cohort %in% "Detroit",]

T3_EN_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_a_dt[,2]>=pr,1,0))), factor((T3_EN_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_a_dt[,2]>pr,1,0))), factor(as.character(T3_EN_a_dt[,1])),positive="1")
  T3_EN_tabout_dt=rbind(T3_EN_tabout_dt,c(pr,aa$byClass[1:4]))
}

T3_EN_RC_st=roc(response=T3_EN_a[grep("Stanford", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Stanford", rownames(T3_EN_a)),2],direction="<")
T3_EN_st_auc <- round(T3_EN_RC_st$auc, digits = 2)
T3_EN_st_auc_txt <- "0.70"
T3_EN_st_ci <- ci(T3_EN_RC_st, of = "auc")
T3_EN_st_ci_txt <- "95% CI: 0.43-0.97"

T3_EN_a_st <- T3_EN_a2[T3_EN_a2$Cohort %in% "Stanford",]

T3_EN_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_a_st[,2]>=pr,1,0))), factor((T3_EN_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_a_st[,2]>pr,1,0))), factor(as.character(T3_EN_a_st[,1])),positive="1")
  T3_EN_tabout_st=rbind(T3_EN_tabout_st,c(pr,aa$byClass[1:4]))
}

T3_EN_RC_osl=roc(response=T3_EN_a[grep("Oslo", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Oslo", rownames(T3_EN_a)),2],direction="<")
T3_EN_osl_auc <- round(T3_EN_RC_osl$auc, digits = 2)
T3_EN_osl_ci <- ci(T3_EN_RC_osl, of = "auc")
T3_EN_osl_ci_txt <- "95% CI: 0.49-0.73"

T3_EN_a_osl <- T3_EN_a2[T3_EN_a2$Cohort %in% "Oslo",]

T3_EN_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_a_osl[,2]>=pr,1,0))), factor((T3_EN_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_a_osl[,2]>pr,1,0))), factor(as.character(T3_EN_a_osl[,1])),positive="1")
  T3_EN_tabout_osl=rbind(T3_EN_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
EN_ROC_plot_T3 <- ggroc(list("Detroit" = T3_EN_RC_dt,
                             "Stanford" = T3_EN_RC_st,
                             "Oslo" = T3_EN_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_EN_dt_auc, ', ', T3_EN_dt_ci_txt), 
                    paste0('AUC = ', T3_EN_st_auc_txt, ', ', T3_EN_st_ci_txt), 
                    paste0('AUC = ', T3_EN_osl_auc, ', ', T3_EN_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T3_EN_a_mod <- as.data.frame(T3_EN_a)
T3_EN_a_mod <- merge(T3_EN_a_mod, samp_T3, by = "row.names")

T3_EN_box <- ggplot(T3_EN_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# Elastic net Time interval 3 -------------------------------------------------------------

load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T4_EN.RData")
T4_EN_freq <- freq
T4_EN_pile <- pile

T4_EN_a=NULL
for (i in 1:length(T4_EN_pile)){
  T4_EN_a=rbind(T4_EN_a,T4_EN_pile[[i]]$outmat)
}


T4_EN_a2 <- T4_EN_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")

T4_EN_RC=roc(response=T4_EN_a[,1],predictor=T4_EN_a[,2],direction="<")
T4_EN_RC_dt=roc(response=T4_EN_a[grep("Detroit", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Detroit", rownames(T4_EN_a)),2],direction="<")
T4_EN_dt_auc <- round(T4_EN_RC_dt$auc, digits = 2)
T4_EN_dt_ci <- ci(T4_EN_RC_dt, of = "auc")
T4_EN_dt_ci_txt <- "95% CI: 0.58-0.75"

T4_EN_a_dt <- T4_EN_a2[T4_EN_a2$Cohort %in% "Detroit",]

T4_EN_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_EN_a_dt[,2]>=pr,1,0))), factor((T4_EN_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_EN_a_dt[,2]>pr,1,0))), factor(as.character(T4_EN_a_dt[,1])),positive="1")
  T4_EN_tabout_dt=rbind(T4_EN_tabout_dt,c(pr,aa$byClass[1:4]))
}

T4_EN_RC_st=roc(response=T4_EN_a[grep("Stanford", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Stanford", rownames(T4_EN_a)),2],direction="<")
T4_EN_st_auc <- round(T4_EN_RC_st$auc, digits = 2)
T4_EN_st_ci <- ci(T4_EN_RC_st, of = "auc")
T4_EN_st_ci_txt <- "95% CI: 0.36-0.99"

T4_EN_a_st <- T4_EN_a2[T4_EN_a2$Cohort %in% "Stanford",]

T4_EN_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_EN_a_st[,2]>=pr,1,0))), factor((T4_EN_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_EN_a_st[,2]>pr,1,0))), factor(as.character(T4_EN_a_st[,1])),positive="1")
  T4_EN_tabout_st=rbind(T4_EN_tabout_st,c(pr,aa$byClass[1:4]))
}

T4_EN_RC_osl=roc(response=T4_EN_a[grep("Oslo", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Oslo", rownames(T4_EN_a)),2],direction="<")
T4_EN_osl_auc <- round(T4_EN_RC_osl$auc, digits = 2)
T4_EN_osl_auc_txt <- "0.70"
T4_EN_osl_ci <- ci(T4_EN_RC_osl, of = "auc")
T4_EN_osl_ci_txt <- "95% CI: 0.58-0.81"

T4_EN_a_osl <- T4_EN_a2[T4_EN_a2$Cohort %in% "Oslo",]

T4_EN_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_EN_a_osl[,2]>=pr,1,0))), factor((T4_EN_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_EN_a_osl[,2]>pr,1,0))), factor(as.character(T4_EN_a_osl[,1])),positive="1")
  T4_EN_tabout_osl=rbind(T4_EN_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
EN_ROC_plot_T4 <- ggroc(list("Detroit" = T4_EN_RC_dt,
                             "Stanford" = T4_EN_RC_st,
                             "Oslo" = T4_EN_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_EN_dt_auc, ', ', T4_EN_dt_ci_txt), 
                    paste0('AUC = ', T4_EN_st_auc, ', ', T4_EN_st_ci_txt), 
                    paste0('AUC = ', T4_EN_osl_auc_txt, ', ', T4_EN_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T4_EN_a_mod <- as.data.frame(T4_EN_a)
T4_EN_a_mod <- merge(T4_EN_a_mod, samp_T4, by = "row.names")

T4_EN_box <- ggplot(T4_EN_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))




# Random forest Time interval 1 -------------------------------------------------------------

#Time interval 1 (time point 2)
load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T2_RF.RData")
T2_RF_freq <- freq
T2_RF_pile <- pile

T2_RF_a=NULL
for (i in 1:length(T2_RF_pile)){
  T2_RF_a=rbind(T2_RF_a,T2_RF_pile[[i]]$outmat)
}

T2_RF_a2 <- T2_RF_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")

T2_RF_RC=roc(response=T2_RF_a[,1],predictor=T2_RF_a[,2],direction="<")
T2_RF_RC_dt=roc(response=T2_RF_a[grep("Detroit", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Detroit", rownames(T2_RF_a)),2],direction="<")
T2_RF_dt_auc <- round(T2_RF_RC_dt$auc, digits = 2)
T2_RF_dt_ci <- ci(T2_RF_RC_dt, of = "auc")
T2_RF_dt_ci_txt <- "95% CI: 0.46-0.67"

T2_RF_a_dt <- T2_RF_a2[T2_RF_a2$Cohort %in% "Detroit",]

T2_RF_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_a_dt[,2]>=pr,1,0))), factor((T2_RF_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_a_dt[,2]>pr,1,0))), factor(as.character(T2_RF_a_dt[,1])),positive="1")
  T2_RF_tabout_dt=rbind(T2_RF_tabout_dt,c(pr,aa$byClass[1:4]))
}


T2_RF_RC_st=roc(response=T2_RF_a[grep("Stanford", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Stanford", rownames(T2_RF_a)),2],direction="<")
T2_RF_st_auc <- round(T2_RF_RC_st$auc, digits = 2)
T2_RF_st_ci <- ci(T2_RF_RC_st, of = "auc")
T2_RF_st_ci_txt <- "95% CI: 0.56-0.94"

T2_RF_a_st <- T2_RF_a2[T2_RF_a2$Cohort %in% "Stanford",]

T2_RF_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_a_st[,2]>=pr,1,0))), factor((T2_RF_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_a_st[,2]>pr,1,0))), factor(as.character(T2_RF_a_st[,1])),positive="1")
  T2_RF_tabout_st=rbind(T2_RF_tabout_st,c(pr,aa$byClass[1:4]))
}


T2_RF_RC_osl=roc(response=T2_RF_a[grep("Oslo", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Oslo", rownames(T2_RF_a)),2],direction="<")
T2_RF_osl_auc <- round(T2_RF_RC_osl$auc, digits = 2)
T2_RF_osl_ci <- ci(T2_RF_RC_osl, of = "auc")
T2_RF_osl_ci_txt <- "95% CI: 0.41-0.66"

T2_RF_a_osl <- T2_RF_a2[T2_RF_a2$Cohort %in% "Oslo",]

T2_RF_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_a_osl[,2]>=pr,1,0))), factor((T2_RF_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_a_osl[,2]>pr,1,0))), factor(as.character(T2_RF_a_osl[,1])),positive="1")
  T2_RF_tabout_osl=rbind(T2_RF_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
RF_ROC_plot_T2 <- ggroc(list("Detroit" = T2_RF_RC_dt,
                             "Stanford" = T2_RF_RC_st,
                             "Oslo" = T2_RF_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_RF_dt_auc, ', ', T2_RF_dt_ci_txt), 
                    paste0('AUC = ', T2_RF_st_auc, ', ', T2_RF_st_ci_txt), 
                    paste0('AUC = ', T2_RF_osl_auc, ', ', T2_RF_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T2_RF_a_mod <- as.data.frame(T2_RF_a)
T2_RF_a_mod <- merge(T2_RF_a_mod, samp_T2, by = "row.names")

T2_RF_box <- ggplot(T2_RF_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



# Random forest Time interval 2-------------------------------------------------------------

load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T3_RF.RData")
T3_RF_freq <- freq
T3_RF_pile <- pile

T3_RF_a=NULL
for (i in 1:length(T3_RF_pile)){
  T3_RF_a=rbind(T3_RF_a,T3_RF_pile[[i]]$outmat)
}

T3_RF_a2 <- T3_RF_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")


T3_RF_RC=roc(response=T3_RF_a[,1],predictor=T3_RF_a[,2],direction="<")
T3_RF_RC_dt=roc(response=T3_RF_a[grep("Detroit", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Detroit", rownames(T3_RF_a)),2],direction="<")
T3_RF_dt_auc <- round(T3_RF_RC_dt$auc, digits = 2)
T3_RF_dt_ci <- ci(T3_RF_RC_dt, of = "auc")
T3_RF_dt_ci_txt <- "95% CI: 0.47-0.62"

T3_RF_a_dt <- T3_RF_a2[T3_RF_a2$Cohort %in% "Detroit",]

T3_RF_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_a_dt[,2]>=pr,1,0))), factor((T3_RF_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_a_dt[,2]>pr,1,0))), factor(as.character(T3_RF_a_dt[,1])),positive="1")
  T3_RF_tabout_dt=rbind(T3_RF_tabout_dt,c(pr,aa$byClass[1:4]))
}

T3_RF_RC_st=roc(response=T3_RF_a[grep("Stanford", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Stanford", rownames(T3_RF_a)),2],direction="<")
T3_RF_st_auc <- round(T3_RF_RC_st$auc, digits = 2)
T3_RF_st_ci <- ci(T3_RF_RC_st, of = "auc")
T3_RF_st_ci_txt <- "95% CI: 0.67-1"

T3_RF_a_st <- T3_RF_a2[T3_RF_a2$Cohort %in% "Stanford",]

T3_RF_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_a_st[,2]>=pr,1,0))), factor((T3_RF_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_a_st[,2]>pr,1,0))), factor(as.character(T3_RF_a_st[,1])),positive="1")
  T3_RF_tabout_st=rbind(T3_RF_tabout_st,c(pr,aa$byClass[1:4]))
}

T3_RF_RC_osl=roc(response=T3_RF_a[grep("Oslo", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Oslo", rownames(T3_RF_a)),2],direction="<")
T3_RF_osl_auc <- round(T3_RF_RC_osl$auc, digits = 2)
T3_RF_osl_ci <- ci(T3_RF_RC_osl, of = "auc")
T3_RF_osl_ci_txt <- "95% CI: 0.35-0.6"

T3_RF_a_osl <- T3_RF_a2[T3_RF_a2$Cohort %in% "Oslo",]

T3_RF_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_a_osl[,2]>=pr,1,0))), factor((T3_RF_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_a_osl[,2]>pr,1,0))), factor(as.character(T3_RF_a_osl[,1])),positive="1")
  T3_RF_tabout_osl=rbind(T3_RF_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
RF_ROC_plot_T3 <- ggroc(list("Detroit" = T3_RF_RC_dt,
                             "Stanford" = T3_RF_RC_st,
                             "Oslo" = T3_RF_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_RF_dt_auc, ', ', T3_RF_dt_ci_txt), 
                    paste0('AUC = ', T3_RF_st_auc, ', ', T3_RF_st_ci_txt), 
                    paste0('AUC = ', T3_RF_osl_auc, ', ', T3_RF_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T3_RF_a_mod <- as.data.frame(T3_RF_a)
T3_RF_a_mod <- merge(T3_RF_a_mod, samp_T3, by = "row.names")

T3_RF_box <- ggplot(T3_RF_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# Random forest Time interval 3 -------------------------------------------------------------

load("../Prediction/Predict_2vs1_all_proteins/Predict_2vs1_all_proteins_T4_RF.RData")
T4_RF_freq <- freq
T4_RF_pile <- pile

T4_RF_a=NULL
for (i in 1:length(T4_RF_pile)){
  T4_RF_a=rbind(T4_RF_a,T4_RF_pile[[i]]$outmat)
}

T4_RF_a2 <- T4_RF_a %>%
  as.data.frame() %>%
  rownames_to_column(var = "rownames") %>%
  mutate("rownames2" = rownames) %>%
  column_to_rownames(var = "rownames") %>%
  separate(rownames2, "Cohort")


T4_RF_RC=roc(response=T4_RF_a[,1],predictor=T4_RF_a[,2],direction="<")
T4_RF_RC_dt=roc(response=T4_RF_a[grep("Detroit", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Detroit", rownames(T4_RF_a)),2],direction="<")
T4_RF_dt_auc <- round(T4_RF_RC_dt$auc, digits = 2)
T4_RF_dt_ci <- ci(T4_RF_RC_dt, of = "auc")
T4_RF_dt_ci_txt <- "95% CI: 0.68-0.83"

T4_RF_a_dt <- T4_RF_a2[T4_RF_a2$Cohort %in% "Detroit",]

T4_RF_tabout_dt=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_RF_a_dt[,2]>=pr,1,0))), factor((T4_RF_a_dt[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_RF_a_dt[,2]>pr,1,0))), factor(as.character(T4_RF_a_dt[,1])),positive="1")
  T4_RF_tabout_dt=rbind(T4_RF_tabout_dt,c(pr,aa$byClass[1:4]))
}

T4_RF_RC_st=roc(response=T4_RF_a[grep("Stanford", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Stanford", rownames(T4_RF_a)),2],direction="<")
T4_RF_st_auc <- round(T4_RF_RC_st$auc, digits = 2)
T4_RF_st_ci <- ci(T4_RF_RC_st, of = "auc")
T4_RF_st_ci_txt <- "95% CI: 0.91-1.00"

T4_RF_a_st <- T4_RF_a2[T4_RF_a2$Cohort %in% "Stanford",]

T4_RF_tabout_st=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_RF_a_st[,2]>=pr,1,0))), factor((T4_RF_a_st[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_RF_a_st[,2]>pr,1,0))), factor(as.character(T4_RF_a_st[,1])),positive="1")
  T4_RF_tabout_st=rbind(T4_RF_tabout_st,c(pr,aa$byClass[1:4]))
}

T4_RF_RC_osl=roc(response=T4_RF_a[grep("Oslo", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Oslo", rownames(T4_RF_a)),2],direction="<")
T4_RF_osl_auc <- round(T4_RF_RC_osl$auc, digits = 2)
T4_RF_osl_ci <- ci(T4_RF_RC_osl, of = "auc")
T4_RF_osl_ci_txt <- "95% CI: 0.63-0.84"

T4_RF_a_osl <- T4_RF_a2[T4_RF_a2$Cohort %in% "Oslo",]

T4_RF_tabout_osl=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_RF_a_osl[,2]>=pr,1,0))), factor((T4_RF_a_osl[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_RF_a_osl[,2]>pr,1,0))), factor(as.character(T4_RF_a_osl[,1])),positive="1")
  T4_RF_tabout_osl=rbind(T4_RF_tabout_osl,c(pr,aa$byClass[1:4]))
}

#Plot
RF_ROC_plot_T4 <- ggroc(list("Detroit" = T4_RF_RC_dt,
                             "Stanford" = T4_RF_RC_st,
                             "Oslo" = T4_RF_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_RF_dt_auc, ', ', T4_RF_dt_ci_txt), 
                    paste0('AUC = ', T4_RF_st_auc, ', ', T4_RF_st_ci_txt), 
                    paste0('AUC = ', T4_RF_osl_auc, ', ', T4_RF_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


#Boxplot
T4_RF_a_mod <- as.data.frame(T4_RF_a)
T4_RF_a_mod <- merge(T4_RF_a_mod, samp_T4, by = "row.names")

T4_RF_box <- ggplot(T4_RF_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))


# Merge plots -------------------------------------------------------------


get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

for_legend_roc <- ggroc(list("Detroit" = T4_RF_RC_dt,
                             "Stanford" = T4_RF_RC_st,
                             "Oslo" = T4_RF_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  theme_bw() +
  theme(legend.title = element_blank())



legend_roc <- get_legend(for_legend_roc)

for_legend_box <- ggplot(T3_RF_a_mod, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot() +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
   theme(legend.title = element_blank())

legend_box <- get_legend(for_legend_box)



EN_ROC_plots <- grid.arrange(EN_ROC_plot_T2, EN_ROC_plot_T3, EN_ROC_plot_T4, legend_roc,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.3))

EN_box_plots <- grid.arrange(T2_EN_box, T3_EN_box, T4_EN_box, legend_box,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.3))


RF_ROC_plots <- grid.arrange(RF_ROC_plot_T2, RF_ROC_plot_T3, RF_ROC_plot_T4, legend_roc,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.3))

RF_box_plots <- grid.arrange(T2_RF_box, T3_RF_box, T4_RF_box, legend_box,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.3))



title1 <- textGrob("Elastic Net", gp = gpar(fontsize = 20))
title2 <- textGrob("Random Forest", gp = gpar(fontsize = 20))

all_plots <- grid.arrange(title1, EN_ROC_plots, EN_box_plots, 
                          title2, RF_ROC_plots, RF_box_plots,
                          ncol = 1, heights = c(0.2, 1, 1, 0.2, 1, 1))

all_plots <- ggarrange(title1, EN_ROC_plots, EN_box_plots, 
                          title2, RF_ROC_plots, RF_box_plots,
                          ncol = 1, heights = c(0.2, 1, 1, 0.2, 1, 1),
                       labels = c(NA, "A)", "B)", NA, "C)", "D)"))

ggsave(filename= "All_roc_box_2vs1_labels.png",
       plot = all_plots,
       device = "png",
       path = "../Prediction/",
       width = 35,
       height = 45,
       units = "cm",
       bg = "white")





#Arrange figure for main figures
EN_ROC_plot_T4_2 <- ggroc(list("Detroit" = T4_EN_RC_dt,
                               "Stanford" = T4_EN_RC_st,
                               "Oslo" = T4_EN_RC_osl),
                          aes = c("colour", "size"),
                          legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Elastic Net") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_EN_dt_auc, ', ', T4_EN_dt_ci_txt), 
                    paste0('AUC = ', T4_EN_st_auc, ', ', T4_EN_st_ci_txt), 
                    paste0('AUC = ', T4_EN_osl_auc_txt, ', ', T4_EN_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))


RF_ROC_plot_T4_2 <- ggroc(list("Detroit" = T4_RF_RC_dt,
                             "Stanford" = T4_RF_RC_st,
                             "Oslo" = T4_RF_RC_osl),
                        aes = c("colour", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D","#619CFF", "#00BA38")) +  
  scale_size_manual(values=c(0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Random Forest") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_RF_dt_auc, ', ', T4_RF_dt_ci_txt), 
                    paste0('AUC = ', T4_RF_st_auc, ', ', T4_RF_st_ci_txt), 
                    paste0('AUC = ', T4_RF_osl_auc, ', ', T4_RF_osl_ci_txt)),
           color = c("#F8766D", "#619CFF", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10))

all_plots_2 <- grid.arrange(EN_ROC_plot_T4_2, RF_ROC_plot_T4_2, legend_roc,
                          ncol = 3, widths = c(1, 1, 0.3))

ggsave(filename= "All_roc_box_2vs1_main.png",
       plot = all_plots_2,
       device = "png",
       path = "../Prediction/",
       width = 25,
       height = 10,
       units = "cm")

