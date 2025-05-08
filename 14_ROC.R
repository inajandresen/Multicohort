rm(list=ls())
library(pROC)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpattern)
library(ggpubr)
library(caret)

load("../Data/Merged_cohorts_wins_Detrended_exprset.RData")
Merged_cohorts_wins_Detrended_exprset$Group <- gsub("PE", "LOPE", Merged_cohorts_wins_Detrended_exprset$Group)

Merged_cohorts_wins_Detrended_exprset_T2 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T2"]
samp_T2 <- pData(Merged_cohorts_wins_Detrended_exprset_T2)
samp_T2 <- samp_T2 %>% dplyr::select(Group) 
Merged_cohorts_wins_Detrended_exprset_T3 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T3"]
samp_T3 <- pData(Merged_cohorts_wins_Detrended_exprset_T3)
samp_T3 <- samp_T3 %>% dplyr::select(Group)
Merged_cohorts_wins_Detrended_exprset_T4 <- Merged_cohorts_wins_Detrended_exprset[, Merged_cohorts_wins_Detrended_exprset$Timepoint %in% "T4"]
samp_T4 <- pData(Merged_cohorts_wins_Detrended_exprset_T4)
samp_T4 <- samp_T4 %>% dplyr::select(Group)


# T2 (time interval 1) Proteins only ----------------------------------------------------
#All proteins EN
load("../Prediction/Prediction_all_proteins/Predict_T2EN.RData")
T2_EN_freq <- freq
T2_EN_pile <- pile

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}


T2_EN_a_df <- as.data.frame(T2_EN_a)
T2_EN_a_df$out_binary <- ifelse(T2_EN_a_df$out > 0.5, 1, 0)
confusionMatrix(as.factor(T2_EN_a_df$out_binary), as.factor(T2_EN_a_df$outcs))
T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_auc=round(ci.auc(T2_EN_RC),3)

T2_EN_auc <- round(T2_EN_RC$auc, digits = 2)
#T2_EN_auc_auc_txt <- "0.74"
T2_EN_ci <- ci(T2_EN_RC, of = "auc")
T2_EN_ci_txt <- "95% CI: 0.68-0.81"

T2_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_a[,2]>=pr,1,0))), factor((T2_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_a[,2]>pr,1,0))), factor(as.character(T2_EN_a[,1])),positive="1")
  T2_EN_tabout=rbind(T2_EN_tabout,c(pr,aa$byClass[1:4]))
}

#All proteins RF
load("../Prediction/Prediction_all_proteins/Predict_T2RF.RData")
T2_RF_freq <- freq
T2_RF_pile <- pile

T2_RF_a=NULL
for (i in 1:length(T2_RF_pile)){
  T2_RF_a=rbind(T2_RF_a,T2_RF_pile[[i]]$outmat)
}
T2_RF_RC=roc(response=T2_RF_a[,1],predictor=T2_RF_a[,2],direction="<")
T2_RF_auc <- round(T2_RF_RC$auc, digits = 2)
T2_RF_auc_auc_txt <- "0.70"
T2_RF_ci <- ci(T2_RF_RC, of = "auc")
T2_RF_ci_txt <- "95% CI: 0.64-0.77"

T2_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_a[,2]>=pr,1,0))), factor((T2_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_a[,2]>pr,1,0))), factor(as.character(T2_RF_a[,1])),positive="1")
  T2_RF_tabout=rbind(T2_RF_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T2EN.RData")
T2_EN_50_freq <- freq
T2_EN_50_pile <- pile

T2_EN_50_a=NULL
for (i in 1:length(T2_EN_50_pile)){
  T2_EN_50_a=rbind(T2_EN_50_a,T2_EN_50_pile[[i]]$outmat)
}
T2_EN_50_RC=roc(response=T2_EN_50_a[,1],predictor=T2_EN_50_a[,2],direction="<")

T2_EN_50_auc <- round(T2_EN_50_RC$auc, digits = 2)
#T2_EN_50_auc_auc_txt <- "0.67"
T2_EN_50_ci <- ci(T2_EN_50_RC, of = "auc")
T2_EN_50_ci_txt <- "95% CI: 0.60-0.74"

T2_EN_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_EN_50_a[,2]>=pr,1,0))), factor((T2_EN_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_EN_50_a[,2]>pr,1,0))), factor(as.character(T2_EN_50_a[,1])),positive="1")
  T2_EN_50_tabout=rbind(T2_EN_50_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins RF
load("../Prediction/Prediction_top50_wilcoxon/Predict_T2RF.RData")
T2_RF_50_freq <- freq
T2_RF_50_pile <- pile

T2_RF_50_a=NULL
for (i in 1:length(T2_RF_50_pile)){
  T2_RF_50_a=rbind(T2_RF_50_a,T2_RF_50_pile[[i]]$outmat)
}
T2_RF_50_RC=roc(response=T2_RF_50_a[,1],predictor=T2_RF_50_a[,2],direction="<")

T2_RF_50_auc <- round(T2_RF_50_RC$auc, digits = 2)
#T2_RF_50_auc_auc_txt <- "0.72"
T2_RF_50_ci <- ci(T2_RF_50_RC, of = "auc")
T2_RF_50_ci_txt <- "95% CI: 0.65-0.78"

T2_RF_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T2_RF_50_a[,2]>=pr,1,0))), factor((T2_RF_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T2_RF_50_a[,2]>pr,1,0))), factor(as.character(T2_RF_50_a[,1])),positive="1")
  T2_RF_50_tabout=rbind(T2_RF_50_tabout,c(pr,aa$byClass[1:4]))
}

#Plot
ROC_plot_T2 <- ggroc(list("Random Forest (all proteins)" = T2_RF_RC,
                          "Random Forest (50 proteins)" = T2_RF_50_RC,
                          "Elastic net (all proteins)" = T2_EN_RC,
                           "Elastic net (50 proteins)" = T2_EN_50_RC),
                aes = c("colour", 
              "linetype", "size"),
              legacy.axes = TRUE) +
  scale_color_manual(values = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2")) +
  scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 1\nWeek 12-19") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T2_RF_auc_auc_txt, ', ', T2_RF_ci_txt), 
                    paste0('AUC = ', T2_RF_50_auc, ', ', T2_RF_50_ci_txt), 
                    paste0('AUC = ', T2_EN_auc, ', ', T2_EN_ci_txt), 
                    paste0('AUC = ', T2_EN_50_auc, ', ', T2_EN_50_ci_txt)),
           color = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T2_EN_a_mod <- as.data.frame(T2_EN_a)
T2_EN_a_mod$Model <- "Elastic net (All proteins)"
T2_EN_a_mod <- merge(T2_EN_a_mod, samp_T2, by = "row.names")

T2_EN_50_a_mod <- as.data.frame(T2_EN_50_a)
T2_EN_50_a_mod$Model <- "Elastic net (50 proteins)"
T2_EN_50_a_mod <- merge(T2_EN_50_a_mod, samp_T2, by = "row.names")

T2_RF_a_mod <- as.data.frame(T2_RF_a)
T2_RF_a_mod$Model <- "Random Forest (All proteins)"
T2_RF_a_mod <- merge(T2_RF_a_mod, samp_T2, by = "row.names")

T2_RF_50_a_mod <- as.data.frame(T2_RF_50_a)
T2_RF_50_a_mod$Model <- "Random Forest (50 proteins)"
T2_RF_50_a_mod <- merge(T2_RF_50_a_mod, samp_T2, by = "row.names")


T2_EN_RF_res <- rbind(T2_EN_a_mod, T2_EN_50_a_mod, T2_RF_a_mod, T2_RF_50_a_mod)

T2_EN_RF_res$Model <- factor(T2_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (50 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (50 proteins)"))


T2_box <- ggplot(T2_EN_RF_res, aes(x = Group, y=out, fill = Model)) +
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("Random Forest (All proteins)" = "none", 
                                 "Elastic net (All proteins)" = "none",
                                 "Random Forest (50 proteins)" = "stripe",
                                 "Elastic net (50 proteins)" = "stripe")) +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "dodgerblue4", 
                              "Elastic net (All proteins)" = "maroon4",
                              "Random Forest (50 proteins)" = "dodgerblue",
                              "Elastic net (50 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



# T3 (time interval 2) Proteins only ----------------------------------------------------
#All proteins EN
load("../Prediction/Prediction_all_proteins/Predict_T3EN.RData")
T3_EN_freq <- freq
T3_EN_pile <- pile

T3_EN_a=NULL
for (i in 1:length(T3_EN_pile)){
  T3_EN_a=rbind(T3_EN_a,T3_EN_pile[[i]]$outmat)
}
T3_EN_RC=roc(response=T3_EN_a[,1],predictor=T3_EN_a[,2],direction="<")
T3_EN_auc=round(ci.auc(T3_EN_RC),3)

T3_EN_auc <- round(T3_EN_RC$auc, digits = 2)
#T3_EN_auc_auc_txt <- "0.74"
T3_EN_ci <- ci(T3_EN_RC, of = "auc")
T3_EN_ci_txt <- "95% CI: 0.71-0.81"

T3_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_a[,2]>=pr,1,0))), factor((T3_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_a[,2]>pr,1,0))), factor(as.character(T3_EN_a[,1])),positive="1")
  T3_EN_tabout=rbind(T3_EN_tabout,c(pr,aa$byClass[1:4]))
}

#All proteins RF
load("../Prediction/Prediction_all_proteins/Predict_T3RF.RData")
T3_RF_freq <- freq
T3_RF_pile <- pile

T3_RF_a=NULL
for (i in 1:length(T3_RF_pile)){
  T3_RF_a=rbind(T3_RF_a,T3_RF_pile[[i]]$outmat)
}
T3_RF_RC=roc(response=T3_RF_a[,1],predictor=T3_RF_a[,2],direction="<")
T3_RF_auc <- round(T3_RF_RC$auc, digits = 2)
#T3_RF_auc_auc_txt <- "0.76"
T3_RF_ci <- ci(T3_RF_RC, of = "auc")
T3_RF_ci_txt <- "95% CI: 0.70-0.81"

T3_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_a[,2]>=pr,1,0))), factor((T3_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_a[,2]>pr,1,0))), factor(as.character(T3_RF_a[,1])),positive="1")
  T3_RF_tabout=rbind(T3_RF_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T3EN.RData")
T3_EN_50_freq <- freq
T3_EN_50_pile <- pile

T3_EN_50_a=NULL
for (i in 1:length(T3_EN_50_pile)){
  T3_EN_50_a=rbind(T3_EN_50_a,T3_EN_50_pile[[i]]$outmat)
}
T3_EN_50_RC=roc(response=T3_EN_50_a[,1],predictor=T3_EN_50_a[,2],direction="<")

T3_EN_50_auc <- round(T3_EN_50_RC$auc, digits = 2)
#T3_EN_50_auc_auc_txt <- "0.72"
T3_EN_50_ci <- ci(T3_EN_50_RC, of = "auc")
T3_EN_50_ci_txt <- "95% CI: 0.66-0.77"

T3_EN_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_EN_50_a[,2]>=pr,1,0))), factor((T3_EN_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_EN_50_a[,2]>pr,1,0))), factor(as.character(T3_EN_50_a[,1])),positive="1")
  T3_EN_50_tabout=rbind(T3_EN_50_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins RF
load("../Prediction/Prediction_top50_wilcoxon/Predict_T3RF.RData")
T3_RF_50_freq <- freq
T3_RF_50_pile <- pile

T3_RF_50_a=NULL
for (i in 1:length(T3_RF_50_pile)){
  T3_RF_50_a=rbind(T3_RF_50_a,T3_RF_50_pile[[i]]$outmat)
}
T3_RF_50_RC=roc(response=T3_RF_50_a[,1],predictor=T3_RF_50_a[,2],direction="<")

T3_RF_50_auc <- round(T3_RF_50_RC$auc, digits = 2)
#T3_RF_50_auc_auc_txt <- "0.76"
T3_RF_50_ci <- ci(T3_RF_50_RC, of = "auc")
T3_RF_50_ci_txt <- "95% CI: 0.71-0.81"

T3_RF_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T3_RF_50_a[,2]>=pr,1,0))), factor((T3_RF_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T3_RF_50_a[,2]>pr,1,0))), factor(as.character(T3_RF_50_a[,1])),positive="1")
  T3_RF_50_tabout=rbind(T3_RF_50_tabout,c(pr,aa$byClass[1:4]))
}

#Plot
ROC_plot_T3 <- ggroc(list("Random Forest (all proteins)" = T3_RF_RC,
                          "Random Forest (50 proteins)" = T3_RF_50_RC,
                          "Elastic net (all proteins)" = T3_EN_RC,
                          "Elastic net (50 proteins)" = T3_EN_50_RC),
                     aes = c("colour", 
                             "linetype", "size"),
                     legacy.axes = TRUE) +
  scale_color_manual(values = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2")) +
  scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 2\nWeek 19-27") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T3_RF_auc, ', ', T3_RF_ci_txt), 
                    paste0('AUC = ', T3_RF_50_auc, ', ', T3_RF_50_ci_txt), 
                    paste0('AUC = ', T3_EN_auc, ', ', T3_EN_ci_txt), 
                    paste0('AUC = ', T3_EN_50_auc, ', ', T3_EN_50_ci_txt)),
           color = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T3_EN_a_mod <- as.data.frame(T3_EN_a)
T3_EN_a_mod$Model <- "Elastic net (All proteins)"
T3_EN_a_mod <- merge(T3_EN_a_mod, samp_T3, by = "row.names")

T3_EN_50_a_mod <- as.data.frame(T3_EN_50_a)
T3_EN_50_a_mod$Model <- "Elastic net (50 proteins)"
T3_EN_50_a_mod <- merge(T3_EN_50_a_mod, samp_T3, by = "row.names")

T3_RF_a_mod <- as.data.frame(T3_RF_a)
T3_RF_a_mod$Model <- "Random Forest (All proteins)"
T3_RF_a_mod <- merge(T3_RF_a_mod, samp_T3, by = "row.names")

T3_RF_50_a_mod <- as.data.frame(T3_RF_50_a)
T3_RF_50_a_mod$Model <- "Random Forest (50 proteins)"
T3_RF_50_a_mod <- merge(T3_RF_50_a_mod, samp_T3, by = "row.names")


T3_EN_RF_res <- rbind(T3_EN_a_mod, T3_EN_50_a_mod, T3_RF_a_mod, T3_RF_50_a_mod)

T3_EN_RF_res$Model <- factor(T3_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (50 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (50 proteins)"))

T3_box <- ggplot(T3_EN_RF_res, aes(x = Group, y=out, fill = Model)) +
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("Random Forest (All proteins)" = "none", 
                                 "Elastic net (All proteins)" = "none",
                                 "Random Forest (50 proteins)" = "stripe",
                                 "Elastic net (50 proteins)" = "stripe")) +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "dodgerblue4", 
                              "Elastic net (All proteins)" = "maroon4",
                              "Random Forest (50 proteins)" = "dodgerblue",
                              "Elastic net (50 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



# T4 (time interval 4) Proteins only ----------------------------------------------------
#All proteins EN
load("../Prediction/Prediction_all_proteins/Predict_T4EN.RData")
T4_EN_freq <- freq
T4_EN_pile <- pile

T4_EN_a=NULL
for (i in 1:length(T4_EN_pile)){
  T4_EN_a=rbind(T4_EN_a,T4_EN_pile[[i]]$outmat)
}
T4_EN_RC=roc(response=T4_EN_a[,1],predictor=T4_EN_a[,2],direction="<")
T4_EN_auc=round(ci.auc(T4_EN_RC),3)

T4_EN_auc <- round(T4_EN_RC$auc, digits = 2)
#T4_EN_auc_auc_txt <- "0.76"
T4_EN_ci <- ci(T4_EN_RC, of = "auc")
T4_EN_ci_txt <- "95% CI: 0.70-0.81"

T4_EN_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_EN_a[,2]>=pr,1,0))), factor((T4_EN_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_EN_a[,2]>pr,1,0))), factor(as.character(T4_EN_a[,1])),positive="1")
  T4_EN_tabout=rbind(T4_EN_tabout,c(pr,aa$byClass[1:4]))
}

#All proteins RF
load("../Prediction/Prediction_all_proteins/Predict_T4RF.RData")
T4_RF_freq <- freq
T4_RF_pile <- pile

T4_RF_a=NULL
for (i in 1:length(T4_RF_pile)){
  T4_RF_a=rbind(T4_RF_a,T4_RF_pile[[i]]$outmat)
}
T4_RF_RC=roc(response=T4_RF_a[,1],predictor=T4_RF_a[,2],direction="<")
T4_RF_auc <- round(T4_RF_RC$auc, digits = 2)
#T4_RF_auc_auc_txt <- "0.79"
T4_RF_ci <- ci(T4_RF_RC, of = "auc")
T4_RF_ci_txt <- "95% CI: 0.74-0.84"

T4_RF_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_RF_a[,2]>=pr,1,0))), factor((T4_RF_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_RF_a[,2]>pr,1,0))), factor(as.character(T4_RF_a[,1])),positive="1")
  T4_RF_tabout=rbind(T4_RF_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T4EN.RData")
T4_EN_50_freq <- freq
T4_EN_50_pile <- pile

T4_EN_50_a=NULL
for (i in 1:length(T4_EN_50_pile)){
  T4_EN_50_a=rbind(T4_EN_50_a,T4_EN_50_pile[[i]]$outmat)
}
T4_EN_50_RC=roc(response=T4_EN_50_a[,1],predictor=T4_EN_50_a[,2],direction="<")

T4_EN_50_auc <- round(T4_EN_50_RC$auc, digits = 2)
#T4_EN_50_auc_auc_txt <- "0.73"
T4_EN_50_ci <- ci(T4_EN_50_RC, of = "auc")
T4_EN_50_ci_txt <- "95% CI: 0.67-0.79"

T4_EN_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_EN_50_a[,2]>=pr,1,0))), factor((T4_EN_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_EN_50_a[,2]>pr,1,0))), factor(as.character(T4_EN_50_a[,1])),positive="1")
  T4_EN_50_tabout=rbind(T4_EN_50_tabout,c(pr,aa$byClass[1:4]))
}

#50 proteins RF
load("../Prediction/Prediction_top50_wilcoxon/Predict_T4RF.RData")
T4_RF_50_freq <- freq
T4_RF_50_pile <- pile

T4_RF_50_a=NULL
for (i in 1:length(T4_RF_50_pile)){
  T4_RF_50_a=rbind(T4_RF_50_a,T4_RF_50_pile[[i]]$outmat)
}
T4_RF_50_RC=roc(response=T4_RF_50_a[,1],predictor=T4_RF_50_a[,2],direction="<")

T4_RF_50_auc <- round(T4_RF_50_RC$auc, digits = 2)
T4_RF_50_auc_txt <- "0.80"
T4_RF_50_ci <- ci(T4_RF_50_RC, of = "auc")
T4_RF_50_ci_txt <- "95% CI: 0.75-0.85"

T4_RF_50_tabout=NULL
for (targetSP in c(0.9,0.8)){
  prs=seq(0.10,0.80,by=0.005)
  specs=NULL
  for(pr in prs){
    aa=confusionMatrix(data = factor((ifelse(T4_RF_50_a[,2]>=pr,1,0))), factor((T4_RF_50_a[,1])),positive="1")
    specs=c(specs,aa$byClass["Specificity"])
  }
  pr=prs[order(abs(specs-targetSP))][1]
  aa=confusionMatrix(data = factor(as.character(ifelse(T4_RF_50_a[,2]>pr,1,0))), factor(as.character(T4_RF_50_a[,1])),positive="1")
  T4_RF_50_tabout=rbind(T4_RF_50_tabout,c(pr,aa$byClass[1:4]))
}


#Plot
ROC_plot_T4 <- ggroc(list("Random Forest (all proteins)" = T4_RF_RC,
                          "Random Forest (50 proteins)" = T4_RF_50_RC,
                          "Elastic net (all proteins)" = T4_EN_RC,
                          "Elastic net (50 proteins)" = T4_EN_50_RC),
                     aes = c("colour", 
                             "linetype", "size"),
                     legacy.axes = TRUE) +
  scale_color_manual(values = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2")) +
  scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 3\nWeek 27-34") +
  #theme_minimal() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_RF_auc, ', ', T4_RF_ci_txt), 
                    paste0('AUC = ', T4_RF_50_auc_txt, ', ', T4_RF_50_ci_txt), 
                    paste0('AUC = ', T4_EN_auc, ', ', T4_EN_ci_txt), 
                    paste0('AUC = ', T4_EN_50_auc, ', ', T4_EN_50_ci_txt)),
           color = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2"), 
           size = 3, 
           x = 0.62, 
           y = c(0.2, 0.15, 0.10, 0.05))


#Boxplot
T4_EN_a_mod <- as.data.frame(T4_EN_a)
T4_EN_a_mod$Model <- "Elastic net (All proteins)"
T4_EN_a_mod <- merge(T4_EN_a_mod, samp_T4, by = "row.names")

T4_EN_50_a_mod <- as.data.frame(T4_EN_50_a)
T4_EN_50_a_mod$Model <- "Elastic net (50 proteins)"
T4_EN_50_a_mod <- merge(T4_EN_50_a_mod, samp_T4, by = "row.names")

T4_RF_a_mod <- as.data.frame(T4_RF_a)
T4_RF_a_mod$Model <- "Random Forest (All proteins)"
T4_RF_a_mod <- merge(T4_RF_a_mod, samp_T4, by = "row.names")

T4_RF_50_a_mod <- as.data.frame(T4_RF_50_a)
T4_RF_50_a_mod$Model <- "Random Forest (50 proteins)"
T4_RF_50_a_mod <- merge(T4_RF_50_a_mod, samp_T4, by = "row.names")

T4_EN_RF_res <- rbind(T4_EN_a_mod, T4_EN_50_a_mod, T4_RF_a_mod, T4_RF_50_a_mod)

T4_EN_RF_res$Model <- factor(T4_EN_RF_res$Model, levels = c("Random Forest (All proteins)",
                                                            "Random Forest (50 proteins)",
                                                            "Elastic net (All proteins)",
                                                            "Elastic net (50 proteins)"))


T4_box <- ggplot(T4_EN_RF_res, aes(x = Group, y=out, fill = Model)) +
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("Random Forest (All proteins)" = "none", 
                                 "Elastic net (All proteins)" = "none",
                                 "Random Forest (50 proteins)" = "stripe",
                                 "Elastic net (50 proteins)" = "stripe")) +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "dodgerblue4", 
                              "Elastic net (All proteins)" = "maroon4",
                              "Random Forest (50 proteins)" = "dodgerblue",
                              "Elastic net (50 proteins)" = "maroon2")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))




# Arrange and save plots --------------------------------------------------
get_legend <-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

for_legend_roc <- ggroc(list("Random Forest (all proteins)" = T4_RF_RC,
                             "Random Forest (50 proteins)" = T4_RF_50_RC,
                             "Elastic net (all proteins)" = T4_EN_RC,
                             "Elastic net (50 proteins)" = T4_EN_50_RC),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("dodgerblue4", "dodgerblue", "maroon4", "maroon2")) +
  scale_linetype_manual(values=c(1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.7, 0.7, 0.7, 0.7)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  theme_minimal() +
  theme(legend.title = element_blank()) 


legend_roc <- get_legend(for_legend_roc)

for_legend_box <- ggplot(T4_EN_RF_res, aes(x = Group, y=out, fill = Model)) +
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("Random Forest (All proteins)" = "none", 
                                 "Elastic net (All proteins)" = "none",
                                 "Random Forest (50 proteins)" = "stripe",
                                 "Elastic net (50 proteins)" = "stripe")) +
  scale_fill_manual(values= c("Random Forest (All proteins)" = "dodgerblue4", 
                              "Elastic net (All proteins)" = "maroon4",
                              "Random Forest (50 proteins)" = "dodgerblue",
                              "Elastic net (50 proteins)" = "maroon2")) +
  theme_classic() +
  theme(legend.title = element_blank())

legend_box <- get_legend(for_legend_box)




all_plots <- grid.arrange(ROC_plot_T2, ROC_plot_T3, ROC_plot_T4, legend_roc,
                          T2_box, T3_box, T4_box, legend_box, 
                          ncol = 4, 
                          widths = c(0.8, 0.8, 0.8, 0.5),
                          #heights = c(1.2, 0.8),
                          heights = c(1.3, 1))

ROCplots <- grid.arrange(ROC_plot_T2, ROC_plot_T3, ROC_plot_T4, legend_roc,
                        ncol = 4, 
                        widths = c(0.8, 0.8, 0.8, 0.5))

BOXplots <- grid.arrange(T2_box, T3_box, T4_box, legend_box,
                         ncol = 4, 
                         widths = c(0.8, 0.8, 0.8, 0.5))


all_plots <- ggarrange(ROCplots, BOXplots,
                          nrow = 2, 
                          heights = c(1.3, 1),
                          labels = c("A)", "B)"))


ggsave(filename= "All_roc_box.png",
       plot = all_plots,
       device = "png",
       path = "../Prediction/",
       width = 30,
       height = 18,
       units = "cm",
       bg = "white")

ggsave(filename= "Figure7_2.png",
       plot = all_plots,
       device = "png",
       path = "../Figures",
       width = 30,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 350)


