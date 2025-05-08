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


# Elastic net -------------------------------------------------------------

#Time interval 1 (time point 2)
load("../Prediction/Prediction_all_proteins/Predict_T2EN.RData")
T2_EN_freq <- freq
T2_EN_pile <- pile

T2_EN_a=NULL
for (i in 1:length(T2_EN_pile)){
  T2_EN_a=rbind(T2_EN_a,T2_EN_pile[[i]]$outmat)
}
T2_EN_RC=roc(response=T2_EN_a[,1],predictor=T2_EN_a[,2],direction="<")
T2_EN_RC_dt=roc(response=T2_EN_a[grep("Detroit", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Detroit", rownames(T2_EN_a)),2],direction="<")
T2_EN_dt_auc <- round(T2_EN_RC_dt$auc, digits = 2)
T2_EN_dt_ci <- ci(T2_EN_RC_dt, of = "auc")
T2_EN_dt_ci_txt <- "95% CI: 0.73-0.89"

T2_EN_RC_st=roc(response=T2_EN_a[grep("Stanford", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Stanford", rownames(T2_EN_a)),2],direction="<")
T2_EN_st_auc <- round(T2_EN_RC_st$auc, digits = 2)
T2_EN_st_ci <- ci(T2_EN_RC_st, of = "auc")
T2_EN_st_ci_txt <- "95% CI: 0.54-0.92"

T2_EN_RC_osl=roc(response=T2_EN_a[grep("Oslo", rownames(T2_EN_a)),1],predictor=T2_EN_a[grep("Oslo", rownames(T2_EN_a)),2],direction="<")
T2_EN_osl_auc <- round(T2_EN_RC_osl$auc, digits = 2)
T2_EN_osl_ci <- ci(T2_EN_RC_osl, of = "auc")
T2_EN_osl_ci_txt <- "95% CI: 0.55-0.78"


#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T2EN.RData")
T2_EN_50_freq <- freq
T2_EN_50_pile <- pile

T2_EN_50_a=NULL
for (i in 1:length(T2_EN_50_pile)){
  T2_EN_50_a=rbind(T2_EN_50_a,T2_EN_50_pile[[i]]$outmat)
}
T2_EN_50_RC=roc(response=T2_EN_50_a[,1],predictor=T2_EN_50_a[,2],direction="<")
T2_EN_50_RC_dt=roc(response=T2_EN_50_a[grep("Detroit", rownames(T2_EN_50_a)),1],predictor=T2_EN_50_a[grep("Detroit", rownames(T2_EN_50_a)),2],direction="<")
T2_EN_50_dt_auc <- round(T2_EN_50_RC_dt$auc, digits = 2)
T2_EN_50_dt_auc_txt <- "0.70"
T2_EN_50_dt_ci <- ci(T2_EN_50_RC_dt, of = "auc")
T2_EN_50_dt_ci_txt <- "95% CI: 0.61-0.80"

T2_EN_50_RC_st=roc(response=T2_EN_50_a[grep("Stanford", rownames(T2_EN_50_a)),1],predictor=T2_EN_50_a[grep("Stanford", rownames(T2_EN_50_a)),2],direction="<")
T2_EN_50_st_auc <- round(T2_EN_50_RC_st$auc, digits = 2)
T2_EN_50_st_ci <- ci(T2_EN_50_RC_st, of = "auc")
T2_EN_50_st_ci_txt <- "95% CI: 0.52-0.92"

T2_EN_50_RC_osl=roc(response=T2_EN_50_a[grep("Oslo", rownames(T2_EN_50_a)),1],predictor=T2_EN_50_a[grep("Oslo", rownames(T2_EN_50_a)),2],direction="<")
T2_EN_50_osl_auc <- round(T2_EN_50_RC_osl$auc, digits = 2)
T2_EN_50_osl_ci <- ci(T2_EN_50_RC_osl, of = "auc")
T2_EN_50_osl_ci_txt <- "95% CI: 0.46-0.72"


#Plot
EN_ROC_plot_T2 <- ggroc(list("Detroit (all proteins)" = T2_EN_RC_dt,
                          "Detroit (50 proteins)" = T2_EN_50_RC_dt,
                          "Stanford (all proteins)" = T2_EN_RC_st,
                          "Stanford (50 proteins)" = T2_EN_50_RC_st,
                          "Oslo (all proteins)" = T2_EN_RC_osl,
                           "Oslo (50 proteins)" = T2_EN_50_RC_osl),
                     aes = c("colour", 
                             "linetype", "size"),
                     legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
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
                    paste0('AUC = ', T2_EN_50_dt_auc_txt, ', ', T2_EN_50_dt_ci_txt), 
                    paste0('AUC = ', T2_EN_st_auc, ', ', T2_EN_st_ci_txt), 
                    paste0('AUC = ', T2_EN_50_st_auc, ', ', T2_EN_50_st_ci_txt),
                    paste0('AUC = ', T2_EN_osl_auc, ', ', T2_EN_osl_ci_txt), 
                    paste0('AUC = ', T2_EN_50_osl_auc, ', ', T2_EN_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T2_EN_a_mod <- as.data.frame(T2_EN_a)
T2_EN_a_mod$Model <- "All proteins"
T2_EN_a_mod <- merge(T2_EN_a_mod, samp_T2, by = "row.names")

T2_EN_50_a_mod <- as.data.frame(T2_EN_50_a)
T2_EN_50_a_mod$Model <- "50 proteins"
T2_EN_50_a_mod <- merge(T2_EN_50_a_mod, samp_T2, by = "row.names")

T2_EN_res <- rbind(T2_EN_a_mod, T2_EN_50_a_mod)

T2_EN_box <- ggplot(T2_EN_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



#Time interval 3 (time point 4)
load("../Prediction/Prediction_all_proteins/Predict_T3EN.RData")
T3_EN_freq <- freq
T3_EN_pile <- pile

T3_EN_a=NULL
for (i in 1:length(T3_EN_pile)){
  T3_EN_a=rbind(T3_EN_a,T3_EN_pile[[i]]$outmat)
}
T3_EN_RC=roc(response=T3_EN_a[,1],predictor=T3_EN_a[,2],direction="<")
T3_EN_RC_dt=roc(response=T3_EN_a[grep("Detroit", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Detroit", rownames(T3_EN_a)),2],direction="<")
T3_EN_dt_auc <- round(T3_EN_RC_dt$auc, digits = 2)
T3_EN_dt_ci <- ci(T3_EN_RC_dt, of = "auc")
T3_EN_dt_ci_txt <- "95% CI: 0.72-0.84"

T3_EN_RC_st=roc(response=T3_EN_a[grep("Stanford", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Stanford", rownames(T3_EN_a)),2],direction="<")
T3_EN_st_auc <- round(T3_EN_RC_st$auc, digits = 2)
T3_EN_st_ci <- ci(T3_EN_RC_st, of = "auc")
T3_EN_st_ci_txt <- "95% CI: 0.55-1"

T3_EN_RC_osl=roc(response=T3_EN_a[grep("Oslo", rownames(T3_EN_a)),1],predictor=T3_EN_a[grep("Oslo", rownames(T3_EN_a)),2],direction="<")
T3_EN_osl_auc <- round(T3_EN_RC_osl$auc, digits = 2)
T3_EN_osl_ci <- ci(T3_EN_RC_osl, of = "auc")
T3_EN_osl_ci_txt <- "95% CI: 0.55-0.77"


#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T3EN.RData")
T3_EN_50_freq <- freq
T3_EN_50_pile <- pile

T3_EN_50_a=NULL
for (i in 1:length(T3_EN_50_pile)){
  T3_EN_50_a=rbind(T3_EN_50_a,T3_EN_50_pile[[i]]$outmat)
}
T3_EN_50_RC=roc(response=T3_EN_50_a[,1],predictor=T3_EN_50_a[,2],direction="<")
T3_EN_50_RC_dt=roc(response=T3_EN_50_a[grep("Detroit", rownames(T3_EN_50_a)),1],predictor=T3_EN_50_a[grep("Detroit", rownames(T3_EN_50_a)),2],direction="<")
T3_EN_50_dt_auc <- round(T3_EN_50_RC_dt$auc, digits = 2)
T3_EN_50_dt_ci <- ci(T3_EN_50_RC_dt, of = "auc")
T3_EN_50_dt_ci_txt <- "95% CI: 0.66-0.80"

T3_EN_50_RC_st=roc(response=T3_EN_50_a[grep("Stanford", rownames(T3_EN_50_a)),1],predictor=T3_EN_50_a[grep("Stanford", rownames(T3_EN_50_a)),2],direction="<")
T3_EN_50_st_auc <- round(T3_EN_50_RC_st$auc, digits = 2)
T3_EN_50_st_ci <- ci(T3_EN_50_RC_st, of = "auc")
T3_EN_50_st_ci_txt <- "95% CI: 0.74-1"

T3_EN_50_RC_osl=roc(response=T3_EN_50_a[grep("Oslo", rownames(T3_EN_50_a)),1],predictor=T3_EN_50_a[grep("Oslo", rownames(T3_EN_50_a)),2],direction="<")
T3_EN_50_osl_auc <- round(T3_EN_50_RC_osl$auc, digits = 2)
T3_EN_50_osl_ci <- ci(T3_EN_50_RC_osl, of = "auc")
T3_EN_50_osl_ci_txt <- "95% CI: 0.49-0.72"


#Plot
EN_ROC_plot_T3 <- ggroc(list("Detroit (all proteins)" = T3_EN_RC_dt,
                             "Detroit (50 proteins)" = T3_EN_50_RC_dt,
                             "Stanford (all proteins)" = T3_EN_RC_st,
                             "Stanford (50 proteins)" = T3_EN_50_RC_st,
                             "Oslo (all proteins)" = T3_EN_RC_osl,
                             "Oslo (50 proteins)" = T3_EN_50_RC_osl),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
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
                    paste0('AUC = ', T3_EN_50_dt_auc, ', ', T3_EN_50_dt_ci_txt), 
                    paste0('AUC = ', T3_EN_st_auc, ', ', T3_EN_st_ci_txt), 
                    paste0('AUC = ', T3_EN_50_st_auc, ', ', T3_EN_50_st_ci_txt),
                    paste0('AUC = ', T3_EN_osl_auc, ', ', T3_EN_osl_ci_txt), 
                    paste0('AUC = ', T3_EN_50_osl_auc, ', ', T3_EN_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))


#Boxplot
T3_EN_a_mod <- as.data.frame(T3_EN_a)
T3_EN_a_mod$Model <- "All proteins"
T3_EN_a_mod <- merge(T3_EN_a_mod, samp_T3, by = "row.names")

T3_EN_50_a_mod <- as.data.frame(T3_EN_50_a)
T3_EN_50_a_mod$Model <- "50 proteins"
T3_EN_50_a_mod <- merge(T3_EN_50_a_mod, samp_T3, by = "row.names")

T3_EN_res <- rbind(T3_EN_a_mod, T3_EN_50_a_mod)

T3_EN_box <- ggplot(T3_EN_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

#Time interval 3 (time point 3)
load("../Prediction/Prediction_all_proteins/Predict_T4EN.RData")
T4_EN_freq <- freq
T4_EN_pile <- pile

T4_EN_a=NULL
for (i in 1:length(T4_EN_pile)){
  T4_EN_a=rbind(T4_EN_a,T4_EN_pile[[i]]$outmat)
}
T4_EN_RC=roc(response=T4_EN_a[,1],predictor=T4_EN_a[,2],direction="<")
T4_EN_RC_dt=roc(response=T4_EN_a[grep("Detroit", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Detroit", rownames(T4_EN_a)),2],direction="<")
T4_EN_dt_auc <- round(T4_EN_RC_dt$auc, digits = 2)
T4_EN_dt_auc_txt <- "0.80"
T4_EN_dt_ci <- ci(T4_EN_RC_dt, of = "auc")
T4_EN_dt_ci_txt <- "95% CI: 0.73-0.87"

T4_EN_RC_st=roc(response=T4_EN_a[grep("Stanford", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Stanford", rownames(T4_EN_a)),2],direction="<")
T4_EN_st_auc <- round(T4_EN_RC_st$auc, digits = 2)
T4_EN_st_auc_txt <- "0.70"
T4_EN_st_ci <- ci(T4_EN_RC_st, of = "auc")
T4_EN_st_ci_txt <- "95% CI: 0.35-1"

T4_EN_RC_osl=roc(response=T4_EN_a[grep("Oslo", rownames(T4_EN_a)),1],predictor=T4_EN_a[grep("Oslo", rownames(T4_EN_a)),2],direction="<")
T4_EN_osl_auc <- round(T4_EN_RC_osl$auc, digits = 2)
T4_EN_osl_auc_txt <- "0.70"
T4_EN_osl_ci <- ci(T4_EN_RC_osl, of = "auc")
T4_EN_osl_ci_txt <- "95% CI: 0.58-0.81"


#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T4EN.RData")
T4_EN_50_freq <- freq
T4_EN_50_pile <- pile

T4_EN_50_a=NULL
for (i in 1:length(T4_EN_50_pile)){
  T4_EN_50_a=rbind(T4_EN_50_a,T4_EN_50_pile[[i]]$outmat)
}
T4_EN_50_RC=roc(response=T4_EN_50_a[,1],predictor=T4_EN_50_a[,2],direction="<")
T4_EN_50_RC_dt=roc(response=T4_EN_50_a[grep("Detroit", rownames(T4_EN_50_a)),1],predictor=T4_EN_50_a[grep("Detroit", rownames(T4_EN_50_a)),2],direction="<")
T4_EN_50_dt_auc <- round(T4_EN_50_RC_dt$auc, digits = 2)
T4_EN_50_dt_ci <- ci(T4_EN_50_RC_dt, of = "auc")
T4_EN_50_dt_ci_txt <- "95% CI: 0.66-0.82"

T4_EN_50_RC_st=roc(response=T4_EN_50_a[grep("Stanford", rownames(T4_EN_50_a)),1],predictor=T4_EN_50_a[grep("Stanford", rownames(T4_EN_50_a)),2],direction="<")
T4_EN_50_st_auc <- round(T4_EN_50_RC_st$auc, digits = 2)
T4_EN_50_st_ci <- ci(T4_EN_50_RC_st, of = "auc")
T4_EN_50_st_ci_txt <- "95% CI: 0.43-1"

T4_EN_50_RC_osl=roc(response=T4_EN_50_a[grep("Oslo", rownames(T4_EN_50_a)),1],predictor=T4_EN_50_a[grep("Oslo", rownames(T4_EN_50_a)),2],direction="<")
T4_EN_50_osl_auc <- round(T4_EN_50_RC_osl$auc, digits = 2)
T4_EN_50_osl_ci <- ci(T4_EN_50_RC_osl, of = "auc")
T4_EN_50_osl_ci_txt <- "95% CI: 0.61-0.83"

#Plot
EN_ROC_plot_T4 <- ggroc(list("Detroit (all proteins)" = T4_EN_RC_dt,
                             "Detroit (50 proteins)" = T4_EN_50_RC_dt,
                             "Stanford (all proteins)" = T4_EN_RC_st,
                             "Stanford (50 proteins)" = T4_EN_50_RC_st,
                             "Oslo (all proteins)" = T4_EN_RC_osl,
                             "Oslo (50 proteins)" = T4_EN_50_RC_osl),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
               color="black", linetype="solid", size = 0.5) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           label= c(paste0('AUC = ', T4_EN_dt_auc_txt, ', ', T4_EN_dt_ci_txt), 
                    paste0('AUC = ', T4_EN_50_dt_auc, ', ', T4_EN_50_dt_ci_txt), 
                    paste0('AUC = ', T4_EN_st_auc_txt, ', ', T4_EN_st_ci_txt), 
                    paste0('AUC = ', T4_EN_50_st_auc, ', ', T4_EN_50_st_ci_txt),
                    paste0('AUC = ', T4_EN_osl_auc_txt, ', ', T4_EN_osl_ci_txt), 
                    paste0('AUC = ', T4_EN_50_osl_auc, ', ', T4_EN_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))



#Boxplot
T4_EN_a_mod <- as.data.frame(T4_EN_a)
T4_EN_a_mod$Model <- "All proteins"
T4_EN_a_mod <- merge(T4_EN_a_mod, samp_T4, by = "row.names")

T4_EN_50_a_mod <- as.data.frame(T4_EN_50_a)
T4_EN_50_a_mod$Model <- "50 proteins"
T4_EN_50_a_mod <- merge(T4_EN_50_a_mod, samp_T4, by = "row.names")

T4_EN_res <- rbind(T4_EN_a_mod, T4_EN_50_a_mod)

T4_EN_box <- ggplot(T4_EN_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 3\nWeek 27-34") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))





# Random forest -----------------------------------------------------------
#Time interval 1 (time point 2)
load("../Prediction/Prediction_all_proteins/Predict_T2RF.RData")
T2_RF_freq <- freq
T2_RF_pile <- pile

T2_RF_a=NULL
for (i in 1:length(T2_RF_pile)){
  T2_RF_a=rbind(T2_RF_a,T2_RF_pile[[i]]$outmat)
}
T2_RF_RC=roc(response=T2_RF_a[,1],predictor=T2_RF_a[,2],direction="<")
T2_RF_RC_dt=roc(response=T2_RF_a[grep("Detroit", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Detroit", rownames(T2_RF_a)),2],direction="<")
T2_RF_dt_auc <- round(T2_RF_RC_dt$auc, digits = 2)
T2_RF_dt_ci <- ci(T2_RF_RC_dt, of = "auc")
T2_RF_dt_ci_txt <- "95% CI: 0.64-0.83"

T2_RF_RC_st=roc(response=T2_RF_a[grep("Stanford", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Stanford", rownames(T2_RF_a)),2],direction="<")
T2_RF_st_auc <- round(T2_RF_RC_st$auc, digits = 2)
T2_RF_st_ci <- ci(T2_RF_RC_st, of = "auc")
T2_RF_st_ci_txt <- "95% CI: 0.56-0.93"

T2_RF_RC_osl=roc(response=T2_RF_a[grep("Oslo", rownames(T2_RF_a)),1],predictor=T2_RF_a[grep("Oslo", rownames(T2_RF_a)),2],direction="<")
T2_RF_osl_auc <- round(T2_RF_RC_osl$auc, digits = 2)
T2_RF_osl_ci <- ci(T2_RF_RC_osl, of = "auc")
T2_RF_osl_ci_txt <- "95% CI: 0.55-0.76"


#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T2RF.RData")
T2_RF_50_freq <- freq
T2_RF_50_pile <- pile

T2_RF_50_a=NULL
for (i in 1:length(T2_RF_50_pile)){
  T2_RF_50_a=rbind(T2_RF_50_a,T2_RF_50_pile[[i]]$outmat)
}
T2_RF_50_RC=roc(response=T2_RF_50_a[,1],predictor=T2_RF_50_a[,2],direction="<")
T2_RF_50_RC_dt=roc(response=T2_RF_50_a[grep("Detroit", rownames(T2_RF_50_a)),1],predictor=T2_RF_50_a[grep("Detroit", rownames(T2_RF_50_a)),2],direction="<")
T2_RF_50_dt_auc <- round(T2_RF_50_RC_dt$auc, digits = 2)
T2_RF_50_dt_ci <- ci(T2_RF_50_RC_dt, of = "auc")
T2_RF_50_dt_ci_txt <- "95% CI: 0.66-0.84"

T2_RF_50_RC_st=roc(response=T2_RF_50_a[grep("Stanford", rownames(T2_RF_50_a)),1],predictor=T2_RF_50_a[grep("Stanford", rownames(T2_RF_50_a)),2],direction="<")
T2_RF_50_st_auc <- round(T2_RF_50_RC_st$auc, digits = 2)
T2_RF_50_st_ci <- ci(T2_RF_50_RC_st, of = "auc")
T2_RF_50_st_ci_txt <- "95% CI: 0.60-0.95"

T2_RF_50_RC_osl=roc(response=T2_RF_50_a[grep("Oslo", rownames(T2_RF_50_a)),1],predictor=T2_RF_50_a[grep("Oslo", rownames(T2_RF_50_a)),2],direction="<")
T2_RF_50_osl_auc <- round(T2_RF_50_RC_osl$auc, digits = 2)
T2_RF_50_osl_ci <- ci(T2_RF_50_RC_osl, of = "auc")
T2_RF_50_osl_ci_txt <- "95% CI: 0.55-0.78"


#Plot
RF_ROC_plot_T2 <- ggroc(list("Detroit (all proteins)" = T2_RF_RC_dt,
                             "Detroit (50 proteins)" = T2_RF_50_RC_dt,
                             "Stanford (all proteins)" = T2_RF_RC_st,
                             "Stanford (50 proteins)" = T2_RF_50_RC_st,
                             "Oslo (all proteins)" = T2_RF_RC_osl,
                             "Oslo (50 proteins)" = T2_RF_50_RC_osl),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
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
                    paste0('AUC = ', T2_RF_50_dt_auc, ', ', T2_RF_50_dt_ci_txt), 
                    paste0('AUC = ', T2_RF_st_auc, ', ', T2_RF_st_ci_txt), 
                    paste0('AUC = ', T2_RF_50_st_auc, ', ', T2_RF_50_st_ci_txt),
                    paste0('AUC = ', T2_RF_osl_auc, ', ', T2_RF_osl_ci_txt), 
                    paste0('AUC = ', T2_RF_50_osl_auc, ', ', T2_RF_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))

#Boxplot
T2_RF_a_mod <- as.data.frame(T2_RF_a)
T2_RF_a_mod$Model <- "All proteins"
T2_RF_a_mod <- merge(T2_RF_a_mod, samp_T2, by = "row.names")

T2_RF_50_a_mod <- as.data.frame(T2_RF_50_a)
T2_RF_50_a_mod$Model <- "50 proteins"
T2_RF_50_a_mod <- merge(T2_RF_50_a_mod, samp_T2, by = "row.names")

T2_RF_res <- rbind(T2_RF_a_mod, T2_RF_50_a_mod)

T2_RF_box <- ggplot(T2_RF_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 1\nWeek 12-19") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))



#Time interval 3 (time point 4)
load("../Prediction/Prediction_all_proteins/Predict_T3RF.RData")
T3_RF_freq <- freq
T3_RF_pile <- pile

T3_RF_a=NULL
for (i in 1:length(T3_RF_pile)){
  T3_RF_a=rbind(T3_RF_a,T3_RF_pile[[i]]$outmat)
}
T3_RF_RC=roc(response=T3_RF_a[,1],predictor=T3_RF_a[,2],direction="<")
T3_RF_RC_dt=roc(response=T3_RF_a[grep("Detroit", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Detroit", rownames(T3_RF_a)),2],direction="<")
T3_RF_dt_auc <- round(T3_RF_RC_dt$auc, digits = 2)
T3_RF_dt_ci <- ci(T3_RF_RC_dt, of = "auc")
T3_RF_dt_ci_txt <- "95% CI: 0.71-0.84"

T3_RF_RC_st=roc(response=T3_RF_a[grep("Stanford", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Stanford", rownames(T3_RF_a)),2],direction="<")
T3_RF_st_auc <- round(T3_RF_RC_st$auc, digits = 2)
T3_RF_st_ci <- ci(T3_RF_RC_st, of = "auc")
T3_RF_st_ci_txt <- "95% CI: 0.68-1"

T3_RF_RC_osl=roc(response=T3_RF_a[grep("Oslo", rownames(T3_RF_a)),1],predictor=T3_RF_a[grep("Oslo", rownames(T3_RF_a)),2],direction="<")
T3_RF_osl_auc <- round(T3_RF_RC_osl$auc, digits = 2)
T3_RF_osl_ci <- ci(T3_RF_RC_osl, of = "auc")
T3_RF_osl_ci_txt <- "95% CI: 0.52-0.75"


#50 proteins EN
load("../Prediction/Prediction_top50_wilcoxon/Predict_T3RF.RData")
T3_RF_50_freq <- freq
T3_RF_50_pile <- pile

T3_RF_50_a=NULL
for (i in 1:length(T3_RF_50_pile)){
  T3_RF_50_a=rbind(T3_RF_50_a,T3_RF_50_pile[[i]]$outmat)
}
T3_RF_50_RC=roc(response=T3_RF_50_a[,1],predictor=T3_RF_50_a[,2],direction="<")
T3_RF_50_RC_dt=roc(response=T3_RF_50_a[grep("Detroit", rownames(T3_RF_50_a)),1],predictor=T3_RF_50_a[grep("Detroit", rownames(T3_RF_50_a)),2],direction="<")
T3_RF_50_dt_auc <- round(T3_RF_50_RC_dt$auc, digits = 2)
T3_RF_50_dt_ci <- ci(T3_RF_50_RC_dt, of = "auc")
T3_RF_50_dt_ci_txt <- "95% CI: 0.73-0.85"

T3_RF_50_RC_st=roc(response=T3_RF_50_a[grep("Stanford", rownames(T3_RF_50_a)),1],predictor=T3_RF_50_a[grep("Stanford", rownames(T3_RF_50_a)),2],direction="<")
T3_RF_50_st_auc <- round(T3_RF_50_RC_st$auc, digits = 2)
T3_RF_50_st_ci <- ci(T3_RF_50_RC_st, of = "auc")
T3_RF_50_st_ci_txt <- "95% CI: 0.59-1"

T3_RF_50_RC_osl=roc(response=T3_RF_50_a[grep("Oslo", rownames(T3_RF_50_a)),1],predictor=T3_RF_50_a[grep("Oslo", rownames(T3_RF_50_a)),2],direction="<")
T3_RF_50_osl_auc <- round(T3_RF_50_RC_osl$auc, digits = 2)
T3_RF_50_osl_ci <- ci(T3_RF_50_RC_osl, of = "auc")
T3_RF_50_osl_ci_txt <- "95% CI: 0.51-0.75"


#Plot
RF_ROC_plot_T3 <- ggroc(list("Detroit (all proteins)" = T3_RF_RC_dt,
                             "Detroit (50 proteins)" = T3_RF_50_RC_dt,
                             "Stanford (all proteins)" = T3_RF_RC_st,
                             "Stanford (50 proteins)" = T3_RF_50_RC_st,
                             "Oslo (all proteins)" = T3_RF_RC_osl,
                             "Oslo (50 proteins)" = T3_RF_50_RC_osl),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
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
                    paste0('AUC = ', T3_RF_50_dt_auc, ', ', T3_RF_50_dt_ci_txt), 
                    paste0('AUC = ', T3_RF_st_auc, ', ', T3_RF_st_ci_txt), 
                    paste0('AUC = ', T3_RF_50_st_auc, ', ', T3_RF_50_st_ci_txt),
                    paste0('AUC = ', T3_RF_osl_auc, ', ', T3_RF_osl_ci_txt), 
                    paste0('AUC = ', T3_RF_50_osl_auc, ', ', T3_RF_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))


#Boxplot
T3_RF_a_mod <- as.data.frame(T3_RF_a)
T3_RF_a_mod$Model <- "All proteins"
T3_RF_a_mod <- merge(T3_RF_a_mod, samp_T3, by = "row.names")

T3_RF_50_a_mod <- as.data.frame(T3_RF_50_a)
T3_RF_50_a_mod$Model <- "50 proteins"
T3_RF_50_a_mod <- merge(T3_RF_50_a_mod, samp_T3, by = "row.names")

T3_RF_res <- rbind(T3_RF_a_mod, T3_RF_50_a_mod)

T3_RF_box <- ggplot(T3_RF_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 2\nWeek 19-27") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

#Time interval 3 (time point 3)
load("../Prediction/Prediction_all_proteins/Predict_T4RF.RData")
T4_RF_freq <- freq
T4_RF_pile <- pile

T4_RF_a=NULL
for (i in 1:length(T4_RF_pile)){
  T4_RF_a=rbind(T4_RF_a,T4_RF_pile[[i]]$outmat)
}
T4_RF_RC=roc(response=T4_RF_a[,1],predictor=T4_RF_a[,2],direction="<")
T4_RF_RC_dt=roc(response=T4_RF_a[grep("Detroit", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Detroit", rownames(T4_RF_a)),2],direction="<")
T4_RF_dt_auc <- round(T4_RF_RC_dt$auc, digits = 2)
T4_RF_dt_ci <- ci(T4_RF_RC_dt, of = "auc")
T4_RF_dt_ci_txt <- "95% CI: 0.68-0.83"

T4_RF_RC_st=roc(response=T4_RF_a[grep("Stanford", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Stanford", rownames(T4_RF_a)),2],direction="<")
T4_RF_st_auc <- round(T4_RF_RC_st$auc, digits = 2)
T4_RF_st_ci <- ci(T4_RF_RC_st, of = "auc")
T4_RF_st_ci_txt <- "95% CI: 0.55-1"

T4_RF_RC_osl=roc(response=T4_RF_a[grep("Oslo", rownames(T4_RF_a)),1],predictor=T4_RF_a[grep("Oslo", rownames(T4_RF_a)),2],direction="<")
T4_RF_osl_auc <- round(T4_RF_RC_osl$auc, digits = 2)
T4_RF_osl_ci <- ci(T4_RF_RC_osl, of = "auc")
T4_RF_osl_ci_txt <- "95% CI: 0.59-0.84"


#50 proteins RF
load("../Prediction/Prediction_top50_wilcoxon/Predict_T4RF.RData")
T4_RF_50_freq <- freq
T4_RF_50_pile <- pile

T4_RF_50_a=NULL
for (i in 1:length(T4_RF_50_pile)){
  T4_RF_50_a=rbind(T4_RF_50_a,T4_RF_50_pile[[i]]$outmat)
}
T4_RF_50_RC=roc(response=T4_RF_50_a[,1],predictor=T4_RF_50_a[,2],direction="<")
T4_RF_50_RC_dt=roc(response=T4_RF_50_a[grep("Detroit", rownames(T4_RF_50_a)),1],predictor=T4_RF_50_a[grep("Detroit", rownames(T4_RF_50_a)),2],direction="<")
T4_RF_50_dt_auc <- round(T4_RF_50_RC_dt$auc, digits = 2)
T4_RF_50_dt_ci <- ci(T4_RF_50_RC_dt, of = "auc")
T4_RF_50_dt_ci_txt <- "95% CI: 0.71-0.86"

T4_RF_50_RC_st=roc(response=T4_RF_50_a[grep("Stanford", rownames(T4_RF_50_a)),1],predictor=T4_RF_50_a[grep("Stanford", rownames(T4_RF_50_a)),2],direction="<")
T4_RF_50_st_auc <- round(T4_RF_50_RC_st$auc, digits = 2)
T4_RF_50_st_ci <- ci(T4_RF_50_RC_st, of = "auc")
T4_RF_50_st_ci_txt <- "95% CI: 0.57-1"

T4_RF_50_RC_osl=roc(response=T4_RF_50_a[grep("Oslo", rownames(T4_RF_50_a)),1],predictor=T4_RF_50_a[grep("Oslo", rownames(T4_RF_50_a)),2],direction="<")
T4_RF_50_osl_auc <- round(T4_RF_50_RC_osl$auc, digits = 2)
T4_RF_50_osl_ci <- ci(T4_RF_50_RC_osl, of = "auc")
T4_RF_50_osl_ci_txt <- "95% CI: 0.62-0.84"

#Plot
RF_ROC_plot_T4 <- ggroc(list("Detroit (all proteins)" = T4_RF_RC_dt,
                             "Detroit (50 proteins)" = T4_RF_50_RC_dt,
                             "Stanford (all proteins)" = T4_RF_RC_st,
                             "Stanford (50 proteins)" = T4_RF_50_RC_st,
                             "Oslo (all proteins)" = T4_RF_RC_osl,
                             "Oslo (50 proteins)" = T4_RF_50_RC_osl),
                        aes = c("colour", 
                                "linetype", "size"),
                        legacy.axes = TRUE) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
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
                    paste0('AUC = ', T4_RF_50_dt_auc, ', ', T4_RF_50_dt_ci_txt), 
                    paste0('AUC = ', T4_RF_st_auc, ', ', T4_RF_st_ci_txt), 
                    paste0('AUC = ', T4_RF_50_st_auc, ', ', T4_RF_50_st_ci_txt),
                    paste0('AUC = ', T4_RF_osl_auc, ', ', T4_RF_osl_ci_txt), 
                    paste0('AUC = ', T4_RF_50_osl_auc, ', ', T4_RF_50_osl_ci_txt)),
           color = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38"), 
           size = 3, 
           x = 0.62, 
           y = c(0.25, 0.2, 0.15, 0.10, 0.05, 0.0))



#Boxplot
T4_RF_a_mod <- as.data.frame(T4_RF_a)
T4_RF_a_mod$Model <- "All proteins"
T4_RF_a_mod <- merge(T4_RF_a_mod, samp_T4, by = "row.names")

T4_RF_50_a_mod <- as.data.frame(T4_RF_50_a)
T4_RF_50_a_mod$Model <- "50 proteins"
T4_RF_50_a_mod <- merge(T4_RF_50_a_mod, samp_T4, by = "row.names")

T4_RF_res <- rbind(T4_RF_a_mod, T4_RF_50_a_mod)

T4_RF_box <- ggplot(T4_RF_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
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

for_legend_roc <- ggroc(list("Detroit (all proteins)" = T4_EN_RC_dt,
                         "Detroit (50 proteins)" = T4_EN_50_RC_dt,
                         "Stanford (all proteins)" = T4_RF_RC_st,
                         "Stanford (50 proteins)" = T4_EN_50_RC_st,
                         "Oslo (all proteins)" = T4_EN_RC_osl,
                         "Oslo (50 proteins)" = T4_RF_50_RC_osl),
                    aes = c("colour", 
                            "linetype", "size")) +
  scale_color_manual(values = c("#F8766D", "#F8766D", "#619CFF", "#619CFF", "#00BA38", "#00BA38")) +  
  scale_linetype_manual(values=c(1, 2, 1, 2, 1, 2), ) +
  scale_size_manual(values=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)) +
  theme_minimal() +
  theme(legend.title = element_blank())


legend_roc <- get_legend(for_legend_roc)

for_legend_box <- ggplot(T4_RF_res, aes(x=Group, y=out, fill = Cohort)) + 
  geom_boxplot_pattern(aes(pattern = Model)) +
  scale_pattern_manual(values= c("All proteins" = "none", 
                                 "50 proteins" = "stripe")) +
  scale_fill_manual(values= c("Detroit" = "#F8766D", 
                              "Oslo" = "#00BA38",
                              "Stanford" = "#619CFF")) +
  theme_classic() +
  xlab("") +
  ylab("Predicted probabilities") +
  ggtitle("Time interval 3") +
  theme(legend.position = "right")

legend_box <- get_legend(for_legend_box)



EN_ROC_plots <- grid.arrange(EN_ROC_plot_T2, EN_ROC_plot_T3, EN_ROC_plot_T4, legend_roc,
                          ncol = 4, widths = c(0.8, 0.8, 0.8, 0.6))

EN_box_plots <- grid.arrange(T2_EN_box, T3_EN_box, T4_EN_box, legend_box,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.6))


RF_ROC_plots <- grid.arrange(RF_ROC_plot_T2, RF_ROC_plot_T3, RF_ROC_plot_T4, legend_roc,
                         ncol = 4, widths = c(0.8, 0.8, 0.8, 0.6))

RF_box_plots <- grid.arrange(T2_RF_box, T3_RF_box, T4_RF_box, legend_box,
                             ncol = 4, widths = c(0.8, 0.8, 0.8, 0.6))



title1 <- textGrob("Elastic Net", gp = gpar(fontsize = 20))
title2 <- textGrob("Random Forest", gp = gpar(fontsize = 20))

all_plots <- grid.arrange(title1, EN_ROC_plots, EN_box_plots, 
                          title2, RF_ROC_plots, RF_box_plots,
                          ncol = 1, heights = c(0.2, 1, 1, 0.2, 1, 1))

all_plots <- ggarrange(title1, EN_ROC_plots, EN_box_plots, 
                          title2, RF_ROC_plots, RF_box_plots,
                          ncol = 1, heights = c(0.2, 1, 1, 0.2, 1, 1),
                       labels = c(NA, "A)", "B)", NA, "C)", "D)"))


ggsave(filename= "all_roc_box_sepcohorts_labels.png",
       plot = all_plots,
       device = "png",
       path = "../Prediction/",
       width = 40,
       height = 45,
       units = "cm",
       bg = "white")

